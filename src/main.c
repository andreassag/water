#define _POSIX_C_SOURCE 200809L

#include "sw.h"

#include <errno.h>
#include <getopt.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef enum {
    MODE_AUTO = 0,
    MODE_SCALAR,
    MODE_SIMD
} RunMode;

typedef SWResult (*AlignFn)(const char *, size_t, const char *, size_t, const SWParams *);

typedef struct {
    long long score_sum;
    int best_score;
    long double cells;
    size_t alignments;
} RunStats;

static void print_usage(const char *prog) {
    fprintf(stderr,
            "Usage: %s --query FILE --target FILE [options]\n"
            "\n"
            "Input files can be FASTA or FASTQ.\n"
            "If files contain multiple records, use --all to score all pairwise combinations.\n"
            "\n"
            "Options:\n"
            "  -q, --query FILE         Query FASTA/FASTQ file\n"
            "  -t, --target FILE        Target FASTA/FASTQ file\n"
            "  -m, --mode MODE          Alignment mode: auto|scalar|simd (default: auto)\n"
            "      --match N            Match score (default: 2)\n"
            "      --mismatch N         Mismatch score (default: -2)\n"
            "      --gap-open N         Affine gap open penalty > 0 (default: 3)\n"
            "      --gap-extend N       Affine gap extend penalty > 0 (default: 1)\n"
            "  -r, --repeat N           Repeat full workload N times for benchmarking (default: 1)\n"
            "  -a, --all                Score all query-target sequence pairs\n"
            "  -A, --show-alignment     Print traceback alignment for the first query-target pair\n"
            "  -h, --help               Show this help\n",
            prog);
}

static int parse_int(const char *arg, int *out) {
    char *end = NULL;
    errno = 0;
    long value = strtol(arg, &end, 10);
    if (errno != 0 || !end || *end != '\0' || value < INT_MIN || value > INT_MAX) {
        return -1;
    }
    *out = (int)value;
    return 0;
}

static int parse_mode(const char *arg, RunMode *mode) {
    if (strcmp(arg, "auto") == 0) {
        *mode = MODE_AUTO;
        return 0;
    }
    if (strcmp(arg, "scalar") == 0) {
        *mode = MODE_SCALAR;
        return 0;
    }
    if (strcmp(arg, "simd") == 0) {
        *mode = MODE_SIMD;
        return 0;
    }
    return -1;
}

static double now_seconds(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1e-9;
}

static RunStats run_once(const SequenceSet *queries,
                         const SequenceSet *targets,
                         int all_pairs,
                         AlignFn fn,
                         const SWParams *params) {
    RunStats stats;
    stats.score_sum = 0;
    stats.best_score = 0;
    stats.cells = 0.0L;
    stats.alignments = 0;

    if (all_pairs) {
        for (size_t qi = 0; qi < queries->count; ++qi) {
            for (size_t ti = 0; ti < targets->count; ++ti) {
                const Sequence *q = &queries->items[qi];
                const Sequence *t = &targets->items[ti];
                SWResult res = fn(q->seq, q->len, t->seq, t->len, params);
                if (res.score > stats.best_score) {
                    stats.best_score = res.score;
                }
                stats.score_sum += res.score;
                stats.cells += (long double)q->len * (long double)t->len;
                stats.alignments += 1;
            }
        }
        return stats;
    }

    const Sequence *q = &queries->items[0];
    const Sequence *t = &targets->items[0];
    SWResult res = fn(q->seq, q->len, t->seq, t->len, params);

    stats.best_score = res.score;
    stats.score_sum = res.score;
    stats.cells = (long double)q->len * (long double)t->len;
    stats.alignments = 1;
    return stats;
}

int main(int argc, char **argv) {
    const char *query_path = NULL;
    const char *target_path = NULL;
    RunMode requested_mode = MODE_AUTO;
    int repeat = 1;
    int all_pairs = 0;
    int show_alignment = 0;

    SWParams params;
    params.match = 2;
    params.mismatch = -2;
    params.gap_open = 3;
    params.gap_extend = 1;

    enum {
        OPT_MATCH = 1000,
        OPT_MISMATCH,
        OPT_GAP_OPEN,
        OPT_GAP_EXTEND
    };

    static const struct option long_opts[] = {
        {"query", required_argument, NULL, 'q'},
        {"target", required_argument, NULL, 't'},
        {"mode", required_argument, NULL, 'm'},
        {"match", required_argument, NULL, OPT_MATCH},
        {"mismatch", required_argument, NULL, OPT_MISMATCH},
        {"gap-open", required_argument, NULL, OPT_GAP_OPEN},
        {"gap-extend", required_argument, NULL, OPT_GAP_EXTEND},
        {"repeat", required_argument, NULL, 'r'},
        {"all", no_argument, NULL, 'a'},
        {"show-alignment", no_argument, NULL, 'A'},
        {"help", no_argument, NULL, 'h'},
        {NULL, 0, NULL, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "q:t:m:r:aAh", long_opts, NULL)) != -1) {
        switch (opt) {
            case 'q':
                query_path = optarg;
                break;
            case 't':
                target_path = optarg;
                break;
            case 'm':
                if (parse_mode(optarg, &requested_mode) != 0) {
                    fprintf(stderr, "Invalid mode '%s'. Use auto|scalar|simd.\n", optarg);
                    return 2;
                }
                break;
            case 'r':
                if (parse_int(optarg, &repeat) != 0 || repeat <= 0) {
                    fprintf(stderr, "Invalid repeat value '%s'.\n", optarg);
                    return 2;
                }
                break;
            case 'a':
                all_pairs = 1;
                break;
            case 'A':
                show_alignment = 1;
                break;
            case OPT_MATCH:
                if (parse_int(optarg, &params.match) != 0) {
                    fprintf(stderr, "Invalid match score '%s'.\n", optarg);
                    return 2;
                }
                break;
            case OPT_MISMATCH:
                if (parse_int(optarg, &params.mismatch) != 0) {
                    fprintf(stderr, "Invalid mismatch score '%s'.\n", optarg);
                    return 2;
                }
                break;
            case OPT_GAP_OPEN:
                if (parse_int(optarg, &params.gap_open) != 0 || params.gap_open <= 0) {
                    fprintf(stderr, "Invalid gap-open value '%s'.\n", optarg);
                    return 2;
                }
                break;
            case OPT_GAP_EXTEND:
                if (parse_int(optarg, &params.gap_extend) != 0 || params.gap_extend <= 0) {
                    fprintf(stderr, "Invalid gap-extend value '%s'.\n", optarg);
                    return 2;
                }
                break;
            case 'h':
                print_usage(argv[0]);
                return 0;
            default:
                print_usage(argv[0]);
                return 2;
        }
    }

    if (!query_path || !target_path) {
        print_usage(argv[0]);
        return 2;
    }

    SequenceSet queries = {0};
    SequenceSet targets = {0};
    char *err_msg = NULL;

    if (load_sequences(query_path, &queries, &err_msg) != 0) {
        fprintf(stderr, "Error loading query file: %s\n", err_msg ? err_msg : "unknown error");
        free(err_msg);
        return 1;
    }

    if (load_sequences(target_path, &targets, &err_msg) != 0) {
        fprintf(stderr, "Error loading target file: %s\n", err_msg ? err_msg : "unknown error");
        free(err_msg);
        free_sequence_set(&queries);
        return 1;
    }

    if (show_alignment && all_pairs) {
        fprintf(stderr, "--show-alignment cannot be combined with --all.\n");
        free_sequence_set(&queries);
        free_sequence_set(&targets);
        return 1;
    }

    int simd_ok = sw_simd_available();
    const char *simd_backend = sw_simd_backend();
    AlignFn fn = sw_scalar_affine;
    const char *selected_mode_name = "scalar";
    const char *requested_mode_name = "auto";
    const char *active_backend = "scalar";

    switch (requested_mode) {
        case MODE_SCALAR:
            requested_mode_name = "scalar";
            fn = sw_scalar_affine;
            selected_mode_name = "scalar";
            active_backend = "scalar";
            break;
        case MODE_SIMD:
            requested_mode_name = "simd";
            if (!simd_ok) {
                fprintf(stderr, "SIMD mode requested but no SIMD backend is available on this system/build.\n");
                free_sequence_set(&queries);
                free_sequence_set(&targets);
                return 1;
            }
            fn = sw_simd_affine_sse2;
            selected_mode_name = "simd";
            active_backend = simd_backend;
            break;
        case MODE_AUTO:
        default:
            requested_mode_name = "auto";
            if (simd_ok) {
                fn = sw_simd_affine_sse2;
                selected_mode_name = "simd";
                active_backend = simd_backend;
            } else {
                fn = sw_scalar_affine;
                selected_mode_name = "scalar";
                active_backend = "scalar";
            }
            break;
    }

    double t0 = now_seconds();
    long long total_score_sum = 0;
    int best_score = 0;
    RunStats first_stats = {0};

    for (int i = 0; i < repeat; ++i) {
        RunStats stats = run_once(&queries, &targets, all_pairs, fn, &params);
        if (i == 0) {
            first_stats = stats;
        }
        if (stats.best_score > best_score) {
            best_score = stats.best_score;
        }
        total_score_sum += stats.score_sum;
    }
    double elapsed = now_seconds() - t0;

    long double total_cells = first_stats.cells * (long double)repeat;
    double avg_ms = (elapsed * 1000.0) / (double)repeat;
    double gcups = (elapsed > 0.0) ? (double)(total_cells / (long double)elapsed / 1e9L) : 0.0;

    printf("requested_mode=%s selected_mode=%s simd_backend=%s active_backend=%s alignments=%zu repeat=%d best_score=%d score_sum=%lld avg_ms=%.6f gcups=%.6f query_count=%zu target_count=%zu\n",
           requested_mode_name,
           selected_mode_name,
           simd_backend,
           active_backend,
           first_stats.alignments,
           repeat,
           best_score,
           total_score_sum,
           avg_ms,
           gcups,
           queries.count,
           targets.count);

    if (show_alignment) {
        const Sequence *q = &queries.items[0];
        const Sequence *t = &targets.items[0];
        SWAlignment alignment;
        if (sw_scalar_affine_trace(q->seq, q->len, t->seq, t->len, &params, &alignment) != 0) {
            fprintf(stderr, "Failed to compute traceback alignment.\n");
            free_sequence_set(&queries);
            free_sequence_set(&targets);
            return 1;
        }

        printf("alignment_query=%s alignment_target=%s alignment_score=%d query_range=%d:%d target_range=%d:%d\n",
               q->name,
               t->name,
               alignment.score,
               alignment.query_start,
               alignment.query_end,
               alignment.target_start,
               alignment.target_end);
        printf("query_aln=%s\n", alignment.query_aln ? alignment.query_aln : "");
        printf("match_aln=%s\n", alignment.match_line ? alignment.match_line : "");
        printf("target_aln=%s\n", alignment.target_aln ? alignment.target_aln : "");

        sw_free_alignment(&alignment);
    }

    free_sequence_set(&queries);
    free_sequence_set(&targets);
    return 0;
}
