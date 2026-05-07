#include "sw.h"

#include <ctype.h>
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
#define SW_X86 1
#else
#define SW_X86 0
#endif

#if SW_X86 && defined(__SSE2__)
#include <emmintrin.h>
#define SW_HAS_SSE2_IMPL 1
#else
#define SW_HAS_SSE2_IMPL 0
#endif

#if SW_X86 && (defined(__GNUC__) || defined(__clang__))
#include <immintrin.h>
#define SW_HAS_AVX2_IMPL 1
#else
#define SW_HAS_AVX2_IMPL 0
#endif

#if SW_HAS_AVX2_IMPL && !defined(__AVX2__) && (defined(__GNUC__) || defined(__clang__))
#define SW_AVX2_TARGET __attribute__((target("avx2")))
#else
#define SW_AVX2_TARGET
#endif

static int16_t clamp_score_i16(int value) {
    if (value > INT16_MAX) {
        return INT16_MAX;
    }
    if (value < INT16_MIN) {
        return INT16_MIN;
    }
    return (int16_t)value;
}

static int16_t clamp_penalty_i16(int value) {
    if (value < 0) {
        value = -value;
    }
    if (value > INT16_MAX) {
        value = INT16_MAX;
    }
    return (int16_t)value;
}

static int sw_env_feature_disabled(const char *name) {
    const char *value = getenv(name);
    if (!value || value[0] == '\0') {
        return 0;
    }
    return (strcmp(value, "0") != 0) ? 1 : 0;
}

static int sw_cpu_supports_sse2(void) {
    if (sw_env_feature_disabled("SW_DISABLE_SSE2")) {
        return 0;
    }
#if SW_HAS_SSE2_IMPL
#if defined(__GNUC__) && SW_X86
    return __builtin_cpu_supports("sse2") ? 1 : 0;
#else
    return 1;
#endif
#else
    return 0;
#endif
}

static int sw_cpu_supports_avx2(void) {
    if (sw_env_feature_disabled("SW_DISABLE_AVX2")) {
        return 0;
    }
#if SW_HAS_AVX2_IMPL
#if defined(__GNUC__) && SW_X86
    return __builtin_cpu_supports("avx2") ? 1 : 0;
#elif defined(__AVX2__)
    return 1;
#else
    return 0;
#endif
#else
    return 0;
#endif
}

#if SW_HAS_SSE2_IMPL
static __m128i *build_query_profile_128(const char *query, size_t query_len,
                                        size_t seg_len, const SWParams *params) {
    const size_t seg_width = 8;
    __m128i *profile = (__m128i *)malloc(256 * seg_len * sizeof(__m128i));
    if (!profile) {
        return NULL;
    }

    for (int c = 0; c < 256; ++c) {
        __m128i *row = profile + ((size_t)c * seg_len);
        unsigned char uc = (unsigned char)c;

        for (size_t i = 0; i < seg_len; ++i) {
            int16_t lanes[8];
            for (size_t lane = 0; lane < seg_width; ++lane) {
                size_t qpos = lane * seg_len + i;
                if (qpos < query_len) {
                    unsigned char qq = (unsigned char)query[qpos];
                    int score = (toupper(qq) == toupper(uc)) ? params->match : params->mismatch;
                    lanes[lane] = clamp_score_i16(score);
                } else {
                    lanes[lane] = 0;
                }
            }

            row[i] = _mm_set_epi16(
                lanes[7], lanes[6], lanes[5], lanes[4],
                lanes[3], lanes[2], lanes[1], lanes[0]);
        }
    }

    return profile;
}

static int extract_max_score_128(__m128i v) {
    int16_t lanes[8];
    _mm_storeu_si128((__m128i *)lanes, v);

    int best = 0;
    for (size_t i = 0; i < 8; ++i) {
        if (lanes[i] > best) {
            best = lanes[i];
        }
    }
    return best;
}

static SWResult sw_simd_affine_sse2_impl(const char *query, size_t query_len,
                                         const char *target, size_t target_len,
                                         const SWParams *params) {
    SWResult result = {0, -1, -1};

    const size_t seg_width = 8;
    const size_t seg_len = (query_len + seg_width - 1) / seg_width;

    __m128i *profile = build_query_profile_128(query, query_len, seg_len, params);
    __m128i *pv_h_store = (__m128i *)calloc(seg_len, sizeof(__m128i));
    __m128i *pv_h_load = (__m128i *)calloc(seg_len, sizeof(__m128i));
    __m128i *pv_e = (__m128i *)malloc(seg_len * sizeof(__m128i));

    if (!profile || !pv_h_store || !pv_h_load || !pv_e) {
        free(profile);
        free(pv_h_store);
        free(pv_h_load);
        free(pv_e);
        return sw_scalar_affine(query, query_len, target, target_len, params);
    }

    const int16_t gap_open = clamp_penalty_i16(params->gap_open);
    const int16_t gap_extend = clamp_penalty_i16(params->gap_extend);

    const int16_t neg_inf = (int16_t)(INT16_MIN / 2);
    const __m128i v_zero = _mm_setzero_si128();
    const __m128i v_gap_open = _mm_set1_epi16(gap_open);
    const __m128i v_gap_extend = _mm_set1_epi16(gap_extend);
    const __m128i v_neg_inf = _mm_set1_epi16(neg_inf);

    for (size_t i = 0; i < seg_len; ++i) {
        pv_e[i] = v_neg_inf;
    }

    __m128i v_max_score = v_zero;

    for (size_t j = 0; j < target_len; ++j) {
        unsigned char tc = (unsigned char)target[j];
        const __m128i *v_profile = profile + ((size_t)toupper(tc) * seg_len);

        __m128i v_f = v_neg_inf;
        __m128i v_h = _mm_slli_si128(pv_h_store[seg_len - 1], 2);

        __m128i *tmp = pv_h_load;
        pv_h_load = pv_h_store;
        pv_h_store = tmp;

        for (size_t i = 0; i < seg_len; ++i) {
            __m128i v_e = pv_e[i];

            v_h = _mm_adds_epi16(v_h, v_profile[i]);
            v_h = _mm_max_epi16(v_h, v_e);
            v_h = _mm_max_epi16(v_h, v_f);
            v_h = _mm_max_epi16(v_h, v_zero);

            pv_h_store[i] = v_h;
            v_max_score = _mm_max_epi16(v_max_score, v_h);

            __m128i v_h_gap = _mm_subs_epi16(v_h, v_gap_open);

            v_e = _mm_subs_epi16(v_e, v_gap_extend);
            v_e = _mm_max_epi16(v_e, v_h_gap);
            pv_e[i] = v_e;

            v_f = _mm_subs_epi16(v_f, v_gap_extend);
            v_f = _mm_max_epi16(v_f, v_h_gap);

            v_h = pv_h_load[i];
        }

        for (size_t k = 0; k < seg_width; ++k) {
            v_f = _mm_slli_si128(v_f, 2);
            int any_update = 0;

            for (size_t i = 0; i < seg_len; ++i) {
                __m128i v_h2 = pv_h_store[i];
                __m128i v_h_old = v_h2;
                v_h2 = _mm_max_epi16(v_h2, v_f);

                __m128i v_changed = _mm_cmpgt_epi16(v_h2, v_h_old);
                if (_mm_movemask_epi8(v_changed)) {
                    any_update = 1;
                }

                pv_h_store[i] = v_h2;
                v_max_score = _mm_max_epi16(v_max_score, v_h2);

                __m128i v_h_gap = _mm_subs_epi16(v_h2, v_gap_open);
                v_f = _mm_subs_epi16(v_f, v_gap_extend);
                v_f = _mm_max_epi16(v_f, v_h_gap);
            }

            if (!any_update) {
                break;
            }
        }
    }

    result.score = extract_max_score_128(v_max_score);

    free(profile);
    free(pv_h_store);
    free(pv_h_load);
    free(pv_e);

    return result;
}
#endif

#if SW_HAS_AVX2_IMPL
static void *sw_aligned_alloc32(size_t count, size_t elem_size, int zero_init) {
    if (count == 0 || elem_size == 0 || count > ((size_t)-1) / elem_size) {
        return NULL;
    }

    size_t bytes = count * elem_size;
    void *ptr = NULL;
    if (posix_memalign(&ptr, 32, bytes) != 0) {
        return NULL;
    }

    if (zero_init) {
        memset(ptr, 0, bytes);
    }

    return ptr;
}

SW_AVX2_TARGET static __m256i *build_query_profile_256(const char *query, size_t query_len,
                                                        size_t seg_len, const SWParams *params,
                                                        const unsigned char used_symbols[256],
                                                        const __m256i *query_valid) {
    const size_t seg_width = 16;
    __m256i *profile = (__m256i *)sw_aligned_alloc32(256 * seg_len, sizeof(__m256i), 1);
    if (!profile) {
        return NULL;
    }

    __m256i *query_chars = (__m256i *)sw_aligned_alloc32(seg_len, sizeof(__m256i), 0);
    if (!query_chars) {
        free(profile);
        free(query_chars);
        return NULL;
    }

    for (size_t i = 0; i < seg_len; ++i) {
        int16_t lanes[16];
        for (size_t lane = 0; lane < seg_width; ++lane) {
            size_t qpos = lane * seg_len + i;
            if (qpos < query_len) {
                lanes[lane] = (int16_t)toupper((unsigned char)query[qpos]);
            } else {
                lanes[lane] = 0;
            }
        }
        query_chars[i] = _mm256_loadu_si256((const __m256i *)lanes);
    }

    const __m256i v_match = _mm256_set1_epi16(clamp_score_i16(params->match));
    const __m256i v_mismatch = _mm256_set1_epi16(clamp_score_i16(params->mismatch));

    for (int c = 0; c < 256; ++c) {
        if (!used_symbols[c]) {
            continue;
        }

        __m256i *row = profile + ((size_t)c * seg_len);
        const __m256i v_char = _mm256_set1_epi16((int16_t)c);
        for (size_t i = 0; i < seg_len; ++i) {
            __m256i cmp = _mm256_cmpeq_epi16(query_chars[i], v_char);
            __m256i score = _mm256_or_si256(_mm256_and_si256(cmp, v_match),
                                            _mm256_andnot_si256(cmp, v_mismatch));
            row[i] = _mm256_and_si256(score, query_valid[i]);
        }
    }

    free(query_chars);
    return profile;
}

SW_AVX2_TARGET static __m256i *build_query_valid_mask_256(size_t query_len, size_t seg_len) {
    const size_t seg_width = 16;
    __m256i *query_valid = (__m256i *)sw_aligned_alloc32(seg_len, sizeof(__m256i), 0);
    if (!query_valid) {
        return NULL;
    }

    for (size_t i = 0; i < seg_len; ++i) {
        int16_t valid[16];
        for (size_t lane = 0; lane < seg_width; ++lane) {
            size_t qpos = lane * seg_len + i;
            valid[lane] = (qpos < query_len) ? (int16_t)-1 : 0;
        }
        query_valid[i] = _mm256_loadu_si256((const __m256i *)valid);
    }

    return query_valid;
}

SW_AVX2_TARGET static int extract_max_score_256(__m256i v) {
    __m128i lo = _mm256_castsi256_si128(v);
    __m128i hi = _mm256_extracti128_si256(v, 1);
    __m128i m = _mm_max_epi16(lo, hi);

    m = _mm_max_epi16(m, _mm_srli_si128(m, 8));
    m = _mm_max_epi16(m, _mm_srli_si128(m, 4));
    m = _mm_max_epi16(m, _mm_srli_si128(m, 2));

    return (int)(int16_t)_mm_extract_epi16(m, 0);
}

SW_AVX2_TARGET static __m256i shift_left_word_256(__m256i v) {
    /*
     * Shift left by one 16-bit lane across the full 256-bit vector.
     * _mm256_slli_si256 shifts within each 128-bit half, so inject the
     * low-half carry (old lane 7) into high-half lane 0.
     */
    __m256i shifted = _mm256_slli_si256(v, 2);
    __m128i low128 = _mm256_castsi256_si128(v);
    __m128i carry = _mm_srli_si128(low128, 14);
    __m256i carry_hi = _mm256_inserti128_si256(_mm256_setzero_si256(), carry, 1);
    return _mm256_or_si256(shifted, carry_hi);
}

SW_AVX2_TARGET static SWResult sw_simd_affine_avx2_impl(const char *query, size_t query_len,
                                                         const char *target, size_t target_len,
                                                         const SWParams *params) {
    SWResult result = {0, -1, -1};

    const size_t seg_width = 16;
    const size_t seg_len = (query_len + seg_width - 1) / seg_width;

    unsigned char used_symbols[256] = {0};
    unsigned char *target_codes = (unsigned char *)malloc(target_len);
    if (!target_codes) {
        return sw_scalar_affine(query, query_len, target, target_len, params);
    }

    for (size_t j = 0; j < target_len; ++j) {
        unsigned char uc = (unsigned char)toupper((unsigned char)target[j]);
        target_codes[j] = uc;
        used_symbols[uc] = 1;
    }

    __m256i *query_valid = build_query_valid_mask_256(query_len, seg_len);
    __m256i *profile = build_query_profile_256(query, query_len, seg_len, params, used_symbols, query_valid);
    __m256i *pv_h_store = (__m256i *)sw_aligned_alloc32(seg_len, sizeof(__m256i), 1);
    __m256i *pv_h_load = (__m256i *)sw_aligned_alloc32(seg_len, sizeof(__m256i), 1);
    __m256i *pv_e = (__m256i *)sw_aligned_alloc32(seg_len, sizeof(__m256i), 0);

    if (!query_valid || !profile || !pv_h_store || !pv_h_load || !pv_e) {
        free(target_codes);
        free(query_valid);
        free(profile);
        free(pv_h_store);
        free(pv_h_load);
        free(pv_e);
        return sw_scalar_affine(query, query_len, target, target_len, params);
    }

    const int16_t gap_open = clamp_penalty_i16(params->gap_open);
    const int16_t gap_extend = clamp_penalty_i16(params->gap_extend);

    const int16_t neg_inf = (int16_t)(INT16_MIN / 2);
    const __m256i v_zero = _mm256_setzero_si256();
    const __m256i v_gap_open = _mm256_set1_epi16(gap_open);
    const __m256i v_gap_extend = _mm256_set1_epi16(gap_extend);
    const __m256i v_neg_inf = _mm256_set1_epi16(neg_inf);

    for (size_t i = 0; i < seg_len; ++i) {
        pv_e[i] = v_neg_inf;
    }

    __m256i v_max_score = v_zero;

    for (size_t j = 0; j < target_len; ++j) {
        const __m256i *v_profile = profile + ((size_t)target_codes[j] * seg_len);

        __m256i v_f = v_neg_inf;
        __m256i v_h = shift_left_word_256(pv_h_store[seg_len - 1]);

        __m256i *tmp = pv_h_load;
        pv_h_load = pv_h_store;
        pv_h_store = tmp;

        for (size_t i = 0; i < seg_len; ++i) {
            __m256i v_e = pv_e[i];
            __m256i v_valid = query_valid[i];

            v_h = _mm256_adds_epi16(v_h, v_profile[i]);
            v_h = _mm256_max_epi16(v_h, v_e);
            v_h = _mm256_max_epi16(v_h, v_f);
            v_h = _mm256_max_epi16(v_h, v_zero);
            v_h = _mm256_and_si256(v_h, v_valid);

            pv_h_store[i] = v_h;
            v_max_score = _mm256_max_epi16(v_max_score, v_h);

            __m256i v_h_gap = _mm256_subs_epi16(v_h, v_gap_open);

            v_e = _mm256_subs_epi16(v_e, v_gap_extend);
            v_e = _mm256_max_epi16(v_e, v_h_gap);
            v_e = _mm256_or_si256(_mm256_and_si256(v_e, v_valid),
                                  _mm256_andnot_si256(v_valid, v_neg_inf));
            pv_e[i] = v_e;

            v_f = _mm256_subs_epi16(v_f, v_gap_extend);
            v_f = _mm256_max_epi16(v_f, v_h_gap);
            v_f = _mm256_or_si256(_mm256_and_si256(v_f, v_valid),
                                  _mm256_andnot_si256(v_valid, v_neg_inf));

            v_h = pv_h_load[i];
        }

        for (size_t k = 0; k < seg_width; ++k) {
            v_f = shift_left_word_256(v_f);
            int any_update = 0;

            for (size_t i = 0; i < seg_len; ++i) {
                __m256i v_valid = query_valid[i];
                __m256i v_h2 = pv_h_store[i];
                __m256i v_h_old = v_h2;
                v_h2 = _mm256_max_epi16(v_h2, v_f);
                v_h2 = _mm256_and_si256(v_h2, v_valid);

                __m256i v_changed = _mm256_cmpgt_epi16(v_h2, v_h_old);
                if (_mm256_movemask_epi8(v_changed)) {
                    any_update = 1;
                }

                pv_h_store[i] = v_h2;
                v_max_score = _mm256_max_epi16(v_max_score, v_h2);

                __m256i v_h_gap = _mm256_subs_epi16(v_h2, v_gap_open);
                v_f = _mm256_subs_epi16(v_f, v_gap_extend);
                v_f = _mm256_max_epi16(v_f, v_h_gap);
                v_f = _mm256_or_si256(_mm256_and_si256(v_f, v_valid),
                                      _mm256_andnot_si256(v_valid, v_neg_inf));
            }

            if (!any_update) {
                break;
            }
        }
    }

    result.score = extract_max_score_256(v_max_score);

    free(target_codes);
    free(query_valid);
    free(profile);
    free(pv_h_store);
    free(pv_h_load);
    free(pv_e);

    return result;
}
#endif

int sw_simd_available(void) {
    return (sw_cpu_supports_avx2() || sw_cpu_supports_sse2()) ? 1 : 0;
}

const char *sw_simd_backend(void) {
    if (sw_cpu_supports_avx2()) {
        return "avx2";
    }
    if (sw_cpu_supports_sse2()) {
        return "sse2";
    }
    return "none";
}

SWResult sw_simd_affine_sse2(const char *query, size_t query_len,
                             const char *target, size_t target_len,
                             const SWParams *params) {
    SWResult result = {0, -1, -1};

    if (!query || !target || !params || query_len == 0 || target_len == 0) {
        return result;
    }

#if SW_HAS_AVX2_IMPL
    if (sw_cpu_supports_avx2()) {
        return sw_simd_affine_avx2_impl(query, query_len, target, target_len, params);
    }
#endif
#if SW_HAS_SSE2_IMPL
    if (sw_cpu_supports_sse2()) {
        return sw_simd_affine_sse2_impl(query, query_len, target, target_len, params);
    }
#endif
    return sw_scalar_affine(query, query_len, target, target_len, params);
}
