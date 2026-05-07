#define _POSIX_C_SOURCE 200809L

#include "sw.h"

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int set_error(char **err_msg, const char *fmt, ...) {
    if (!err_msg) {
        return -1;
    }

    va_list ap;
    va_start(ap, fmt);
    int needed = vsnprintf(NULL, 0, fmt, ap);
    va_end(ap);
    if (needed < 0) {
        return -1;
    }

    char *buf = (char *)malloc((size_t)needed + 1);
    if (!buf) {
        return -1;
    }

    va_start(ap, fmt);
    vsnprintf(buf, (size_t)needed + 1, fmt, ap);
    va_end(ap);

    *err_msg = buf;
    return -1;
}

static char *dup_or_default(const char *s, const char *fallback) {
    const char *src = s && *s ? s : fallback;
    size_t n = strlen(src);
    char *copy = (char *)malloc(n + 1);
    if (!copy) {
        return NULL;
    }
    memcpy(copy, src, n + 1);
    return copy;
}

static int push_sequence(SequenceSet *set, const char *name,
                         const char *seq, size_t len,
                         char **err_msg) {
    if (!set || !seq || len == 0) {
        return set_error(err_msg, "Invalid or empty sequence record");
    }

    Sequence *grown = (Sequence *)realloc(set->items, (set->count + 1) * sizeof(Sequence));
    if (!grown) {
        return set_error(err_msg, "Out of memory while storing sequences");
    }
    set->items = grown;

    char *name_copy = dup_or_default(name, "sequence");
    if (!name_copy) {
        return set_error(err_msg, "Out of memory while copying sequence name");
    }

    char *seq_copy = (char *)malloc(len + 1);
    if (!seq_copy) {
        free(name_copy);
        return set_error(err_msg, "Out of memory while copying sequence");
    }

    memcpy(seq_copy, seq, len);
    seq_copy[len] = '\0';

    set->items[set->count].name = name_copy;
    set->items[set->count].seq = seq_copy;
    set->items[set->count].len = len;
    set->count += 1;
    return 0;
}

static void chomp(char *line, ssize_t *n) {
    while (*n > 0 && (line[*n - 1] == '\n' || line[*n - 1] == '\r')) {
        line[*n - 1] = '\0';
        *n -= 1;
    }
}

static int append_clean_sequence(char **buf, size_t *len, size_t *cap,
                                 const char *line) {
    for (size_t i = 0; line[i] != '\0'; ++i) {
        unsigned char ch = (unsigned char)line[i];
        if (isspace(ch)) {
            continue;
        }

        if (*len + 1 >= *cap) {
            size_t new_cap = (*cap == 0) ? 128 : (*cap * 2);
            char *grown = (char *)realloc(*buf, new_cap);
            if (!grown) {
                return -1;
            }
            *buf = grown;
            *cap = new_cap;
        }

        (*buf)[*len] = (char)toupper(ch);
        *len += 1;
    }

    if (*buf) {
        (*buf)[*len] = '\0';
    }

    return 0;
}

static int parse_fasta(FILE *fp, SequenceSet *out, char **err_msg) {
    char *line = NULL;
    size_t line_cap = 0;
    ssize_t line_len;
    size_t line_no = 0;

    char *curr_name = NULL;
    char *seq_buf = NULL;
    size_t seq_len = 0;
    size_t seq_cap = 0;

    while ((line_len = getline(&line, &line_cap, fp)) != -1) {
        line_no += 1;
        chomp(line, &line_len);

        if (line_len == 0) {
            continue;
        }

        if (line[0] == '>') {
            if (curr_name) {
                if (push_sequence(out, curr_name, seq_buf ? seq_buf : "", seq_len, err_msg) != 0) {
                    free(curr_name);
                    free(seq_buf);
                    free(line);
                    return -1;
                }
                seq_len = 0;
                if (seq_buf) {
                    seq_buf[0] = '\0';
                }
            }

            free(curr_name);
            curr_name = dup_or_default(line + 1, "sequence");
            if (!curr_name) {
                free(seq_buf);
                free(line);
                return set_error(err_msg, "Out of memory while reading FASTA header");
            }
            continue;
        }

        if (!curr_name) {
            free(seq_buf);
            free(line);
            return set_error(err_msg, "Invalid FASTA in input: missing header before sequence at line %zu", line_no);
        }

        if (append_clean_sequence(&seq_buf, &seq_len, &seq_cap, line) != 0) {
            free(curr_name);
            free(seq_buf);
            free(line);
            return set_error(err_msg, "Out of memory while reading FASTA sequence at line %zu", line_no);
        }
    }

    if (curr_name) {
        if (push_sequence(out, curr_name, seq_buf ? seq_buf : "", seq_len, err_msg) != 0) {
            free(curr_name);
            free(seq_buf);
            free(line);
            return -1;
        }
    }

    free(curr_name);
    free(seq_buf);
    free(line);

    if (out->count == 0) {
        return set_error(err_msg, "No FASTA sequences found in file");
    }

    return 0;
}

static int parse_fastq(FILE *fp, SequenceSet *out, char **err_msg) {
    char *hdr = NULL;
    char *seq = NULL;
    char *plus = NULL;
    char *qual = NULL;
    size_t cap_hdr = 0;
    size_t cap_seq = 0;
    size_t cap_plus = 0;
    size_t cap_qual = 0;
    ssize_t n_hdr;
    size_t line_no = 0;

    while ((n_hdr = getline(&hdr, &cap_hdr, fp)) != -1) {
        line_no += 1;
        chomp(hdr, &n_hdr);
        if (n_hdr == 0) {
            continue;
        }
        if (hdr[0] != '@') {
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return set_error(err_msg, "Invalid FASTQ header at line %zu", line_no);
        }

        ssize_t n_seq = getline(&seq, &cap_seq, fp);
        if (n_seq == -1) {
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return set_error(err_msg, "Truncated FASTQ record after line %zu", line_no);
        }
        line_no += 1;
        chomp(seq, &n_seq);

        ssize_t n_plus = getline(&plus, &cap_plus, fp);
        if (n_plus == -1) {
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return set_error(err_msg, "Truncated FASTQ record after line %zu", line_no);
        }
        line_no += 1;
        chomp(plus, &n_plus);
        if (n_plus == 0 || plus[0] != '+') {
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return set_error(err_msg, "Invalid FASTQ '+' line at line %zu", line_no);
        }

        ssize_t n_qual = getline(&qual, &cap_qual, fp);
        if (n_qual == -1) {
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return set_error(err_msg, "Truncated FASTQ quality at line %zu", line_no);
        }
        line_no += 1;
        chomp(qual, &n_qual);

        char *clean = NULL;
        size_t clean_len = 0;
        size_t clean_cap = 0;
        if (append_clean_sequence(&clean, &clean_len, &clean_cap, seq) != 0) {
            free(clean);
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return set_error(err_msg, "Out of memory while reading FASTQ sequence at line %zu", line_no);
        }

        if (push_sequence(out, hdr + 1, clean ? clean : "", clean_len, err_msg) != 0) {
            free(clean);
            free(hdr);
            free(seq);
            free(plus);
            free(qual);
            return -1;
        }

        free(clean);
    }

    free(hdr);
    free(seq);
    free(plus);
    free(qual);

    if (out->count == 0) {
        return set_error(err_msg, "No FASTQ reads found in file");
    }

    return 0;
}

int load_sequences(const char *path, SequenceSet *out, char **err_msg) {
    if (!path || !out) {
        return set_error(err_msg, "Invalid argument to load_sequences");
    }

    out->items = NULL;
    out->count = 0;

    FILE *fp = fopen(path, "r");
    if (!fp) {
        return set_error(err_msg, "Failed to open '%s'", path);
    }

    int first = EOF;
    while ((first = fgetc(fp)) != EOF) {
        if (!isspace(first)) {
            break;
        }
    }

    if (first == EOF) {
        fclose(fp);
        return set_error(err_msg, "Input file '%s' is empty", path);
    }

    rewind(fp);

    int rc;
    if (first == '>') {
        rc = parse_fasta(fp, out, err_msg);
    } else if (first == '@') {
        rc = parse_fastq(fp, out, err_msg);
    } else {
        rc = set_error(err_msg, "Unsupported format in '%s': expected FASTA or FASTQ", path);
    }

    fclose(fp);
    if (rc != 0) {
        free_sequence_set(out);
    }
    return rc;
}

void free_sequence_set(SequenceSet *set) {
    if (!set || !set->items) {
        if (set) {
            set->count = 0;
        }
        return;
    }

    for (size_t i = 0; i < set->count; ++i) {
        free(set->items[i].name);
        free(set->items[i].seq);
    }

    free(set->items);
    set->items = NULL;
    set->count = 0;
}
