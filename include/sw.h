#ifndef SW_H
#define SW_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    int match;
    int mismatch;
    int gap_open;
    int gap_extend;
} SWParams;

typedef struct {
    int score;
    int end_query;
    int end_target;
} SWResult;

typedef struct {
    int score;
    int query_start;
    int query_end;
    int target_start;
    int target_end;
    char *query_aln;
    char *target_aln;
    char *match_line;
} SWAlignment;

typedef struct {
    char *name;
    char *seq;
    size_t len;
} Sequence;

typedef struct {
    Sequence *items;
    size_t count;
} SequenceSet;

int sw_score_char(char a, char b, const SWParams *params);

SWResult sw_scalar_affine(const char *query, size_t query_len,
                          const char *target, size_t target_len,
                          const SWParams *params);

int sw_scalar_affine_trace(const char *query, size_t query_len,
                           const char *target, size_t target_len,
                           const SWParams *params,
                           SWAlignment *out_alignment);

void sw_free_alignment(SWAlignment *alignment);

int sw_simd_available(void);
const char *sw_simd_backend(void);

SWResult sw_simd_affine_sse2(const char *query, size_t query_len,
                             const char *target, size_t target_len,
                             const SWParams *params);

int load_sequences(const char *path, SequenceSet *out, char **err_msg);
void free_sequence_set(SequenceSet *set);

#endif
