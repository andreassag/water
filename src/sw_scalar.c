#include "sw.h"

#include <ctype.h>
#include <stdlib.h>

static inline int max2(int a, int b) {
    return (a > b) ? a : b;
}

int sw_score_char(char a, char b, const SWParams *params) {
    if (toupper((unsigned char)a) == toupper((unsigned char)b)) {
        return params->match;
    }
    return params->mismatch;
}

SWResult sw_scalar_affine(const char *query, size_t query_len,
                          const char *target, size_t target_len,
                          const SWParams *params) {
    SWResult result = {0, -1, -1};
    const int neg_inf = -1000000000;

    if (!query || !target || !params || query_len == 0 || target_len == 0) {
        return result;
    }

    int *h_prev = (int *)calloc(target_len + 1, sizeof(int));
    int *h_curr = (int *)calloc(target_len + 1, sizeof(int));
    int *f = (int *)malloc((target_len + 1) * sizeof(int));

    if (!h_prev || !h_curr || !f) {
        free(h_prev);
        free(h_curr);
        free(f);
        return result;
    }

    for (size_t j = 0; j <= target_len; ++j) {
        f[j] = neg_inf;
    }

    for (size_t i = 1; i <= query_len; ++i) {
        int e = neg_inf;
        h_curr[0] = 0;
        f[0] = neg_inf;

        for (size_t j = 1; j <= target_len; ++j) {
            int e_open = h_curr[j - 1] - params->gap_open;
            int e_extend = e - params->gap_extend;
            e = max2(e_open, e_extend);

            int f_open = h_prev[j] - params->gap_open;
            int f_extend = f[j] - params->gap_extend;
            f[j] = max2(f_open, f_extend);

            int diag = h_prev[j - 1] + sw_score_char(query[i - 1], target[j - 1], params);
            int h = diag;
            h = max2(h, e);
            h = max2(h, f[j]);
            h = max2(h, 0);
            h_curr[j] = h;

            if (h > result.score) {
                result.score = h;
                result.end_query = (int)i - 1;
                result.end_target = (int)j - 1;
            }
        }

        int *tmp = h_prev;
        h_prev = h_curr;
        h_curr = tmp;
    }

    free(h_prev);
    free(h_curr);
    free(f);

    return result;
}

static inline size_t sw_idx(size_t i, size_t j, size_t cols) {
    return i * cols + j;
}

void sw_free_alignment(SWAlignment *alignment) {
    if (!alignment) {
        return;
    }

    free(alignment->query_aln);
    free(alignment->target_aln);
    free(alignment->match_line);

    alignment->query_aln = NULL;
    alignment->target_aln = NULL;
    alignment->match_line = NULL;
    alignment->score = 0;
    alignment->query_start = -1;
    alignment->query_end = -1;
    alignment->target_start = -1;
    alignment->target_end = -1;
}

int sw_scalar_affine_trace(const char *query, size_t query_len,
                           const char *target, size_t target_len,
                           const SWParams *params,
                           SWAlignment *out_alignment) {
    const int neg_inf = -1000000000;

    if (!query || !target || !params || !out_alignment) {
        return -1;
    }

    out_alignment->score = 0;
    out_alignment->query_start = -1;
    out_alignment->query_end = -1;
    out_alignment->target_start = -1;
    out_alignment->target_end = -1;
    out_alignment->query_aln = NULL;
    out_alignment->target_aln = NULL;
    out_alignment->match_line = NULL;

    if (query_len == 0 || target_len == 0) {
        return -1;
    }

    size_t rows = query_len + 1;
    size_t cols = target_len + 1;
    if (rows > SIZE_MAX / cols) {
        return -1;
    }
    size_t cell_count = rows * cols;

    int *h = (int *)calloc(cell_count, sizeof(int));
    int *e = (int *)malloc(cell_count * sizeof(int));
    int *f = (int *)malloc(cell_count * sizeof(int));
    if (!h || !e || !f) {
        free(h);
        free(e);
        free(f);
        return -1;
    }

    for (size_t idx = 0; idx < cell_count; ++idx) {
        e[idx] = neg_inf;
        f[idx] = neg_inf;
    }

    int best_score = 0;
    size_t best_i = 0;
    size_t best_j = 0;

    for (size_t i = 1; i <= query_len; ++i) {
        for (size_t j = 1; j <= target_len; ++j) {
            size_t idx = sw_idx(i, j, cols);
            size_t idx_left = sw_idx(i, j - 1, cols);
            size_t idx_up = sw_idx(i - 1, j, cols);
            size_t idx_diag = sw_idx(i - 1, j - 1, cols);

            int e_open = h[idx_left] - params->gap_open;
            int e_extend = e[idx_left] - params->gap_extend;
            e[idx] = max2(e_open, e_extend);

            int f_open = h[idx_up] - params->gap_open;
            int f_extend = f[idx_up] - params->gap_extend;
            f[idx] = max2(f_open, f_extend);

            int diag = h[idx_diag] + sw_score_char(query[i - 1], target[j - 1], params);
            int hv = diag;
            hv = max2(hv, e[idx]);
            hv = max2(hv, f[idx]);
            hv = max2(hv, 0);
            h[idx] = hv;

            if (hv > best_score) {
                best_score = hv;
                best_i = i;
                best_j = j;
            }
        }
    }

    size_t trace_cap = query_len + target_len + 1;
    char *trace_q = (char *)malloc(trace_cap);
    char *trace_t = (char *)malloc(trace_cap);
    if (!trace_q || !trace_t) {
        free(trace_q);
        free(trace_t);
        free(h);
        free(e);
        free(f);
        return -1;
    }

    enum {
        TRACE_H = 0,
        TRACE_E = 1,
        TRACE_F = 2
    };

    size_t i = best_i;
    size_t j = best_j;
    size_t trace_len = 0;
    int state = TRACE_H;

    while (i > 0 && j > 0) {
        if (trace_len + 1 >= trace_cap) {
            free(trace_q);
            free(trace_t);
            free(h);
            free(e);
            free(f);
            return -1;
        }

        size_t idx = sw_idx(i, j, cols);

        if (state == TRACE_H) {
            int hv = h[idx];
            if (hv <= 0) {
                break;
            }

            int diag = h[sw_idx(i - 1, j - 1, cols)] + sw_score_char(query[i - 1], target[j - 1], params);
            if (hv == diag) {
                trace_q[trace_len] = query[i - 1];
                trace_t[trace_len] = target[j - 1];
                trace_len += 1;
                i -= 1;
                j -= 1;
                continue;
            }

            if (hv == e[idx]) {
                state = TRACE_E;
                continue;
            }

            if (hv == f[idx]) {
                state = TRACE_F;
                continue;
            }

            trace_q[trace_len] = query[i - 1];
            trace_t[trace_len] = target[j - 1];
            trace_len += 1;
            i -= 1;
            j -= 1;
            continue;
        }

        if (state == TRACE_E) {
            int ev = e[idx];
            int e_open = h[sw_idx(i, j - 1, cols)] - params->gap_open;
            int e_extend = e[sw_idx(i, j - 1, cols)] - params->gap_extend;

            trace_q[trace_len] = '-';
            trace_t[trace_len] = target[j - 1];
            trace_len += 1;
            j -= 1;

            state = (ev == e_open) ? TRACE_H : TRACE_E;
            if (ev != e_open && ev != e_extend) {
                state = TRACE_H;
            }
            continue;
        }

        {
            int fv = f[idx];
            int f_open = h[sw_idx(i - 1, j, cols)] - params->gap_open;
            int f_extend = f[sw_idx(i - 1, j, cols)] - params->gap_extend;

            trace_q[trace_len] = query[i - 1];
            trace_t[trace_len] = '-';
            trace_len += 1;
            i -= 1;

            state = (fv == f_open) ? TRACE_H : TRACE_F;
            if (fv != f_open && fv != f_extend) {
                state = TRACE_H;
            }
        }
    }

    char *aligned_q = (char *)malloc(trace_len + 1);
    char *aligned_t = (char *)malloc(trace_len + 1);
    char *match_line = (char *)malloc(trace_len + 1);
    if (!aligned_q || !aligned_t || !match_line) {
        free(aligned_q);
        free(aligned_t);
        free(match_line);
        free(trace_q);
        free(trace_t);
        free(h);
        free(e);
        free(f);
        return -1;
    }

    for (size_t k = 0; k < trace_len; ++k) {
        char qc = trace_q[trace_len - 1 - k];
        char tc = trace_t[trace_len - 1 - k];
        aligned_q[k] = qc;
        aligned_t[k] = tc;

        if (qc == tc && qc != '-') {
            match_line[k] = '|';
        } else if (qc == '-' || tc == '-') {
            match_line[k] = ' ';
        } else {
            match_line[k] = '.';
        }
    }

    aligned_q[trace_len] = '\0';
    aligned_t[trace_len] = '\0';
    match_line[trace_len] = '\0';

    out_alignment->score = best_score;
    out_alignment->query_start = (best_score > 0) ? (int)i : -1;
    out_alignment->query_end = (best_score > 0) ? (int)best_i - 1 : -1;
    out_alignment->target_start = (best_score > 0) ? (int)j : -1;
    out_alignment->target_end = (best_score > 0) ? (int)best_j - 1 : -1;
    out_alignment->query_aln = aligned_q;
    out_alignment->target_aln = aligned_t;
    out_alignment->match_line = match_line;

    free(trace_q);
    free(trace_t);
    free(h);
    free(e);
    free(f);
    return 0;
}
