# SIMD Smith-Waterman in C

This project implements local sequence alignment (Smith-Waterman) with affine gap penalties in two versions:

- Scalar dynamic programming baseline.
- SIMD-optimized striped implementation with runtime-dispatched AVX2 or SSE2 intrinsics.

Both versions support the same scoring parameters and can be selected from the CLI.

## Features

- Affine gap model:
  - gap-open penalty
  - gap-extend penalty
- Input from FASTA or FASTQ files.
- CLI mode selection:
  - `scalar`
  - `simd`
  - `auto` (SIMD when available, otherwise scalar)
- Benchmark script to compare scalar vs SIMD throughput.
- Runtime SIMD backend reporting (`simd_backend`, `active_backend`) in CLI output.

## Build

```bash
make
```

Binary output:

- `build/water`

## CLI Usage

```bash
./build/water --query <query.fa|fq> --target <target.fa|fq> [options]
```

Options:

- `--mode auto|scalar|simd` (default: `auto`)
- `--match <int>` (default: `2`)
- `--mismatch <int>` (default: `-2`)
- `--gap-open <int>` (default: `3`)
- `--gap-extend <int>` (default: `1`)
- `--repeat <N>` for benchmarking loops
- `--all` to align all pairs between multi-record query/target files
- `--show-alignment` to print traceback alignment strings for the first query-target pair

Example:

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode auto
```

To print the aligned sequences:

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode auto --show-alignment
```

Output is a single key=value line suitable for scripts, for example:

```text
requested_mode=auto selected_mode=simd simd_backend=avx2 active_backend=avx2 alignments=1 repeat=1 best_score=10 score_sum=10 avg_ms=0.012345 gcups=0.005184 query_count=1 target_count=1
alignment_query=query1 alignment_target=target1 alignment_score=10 query_range=0:7 target_range=0:7
query_aln=ACACACTA
match_aln=||.|||||
target_aln=AGCACACA
```

## Benchmark Script

Run scalar vs SIMD on generated random DNA sequences:

```bash
./scripts/benchmark_compare.sh
```

Environment overrides:

- `REPEAT` (default `200`)
- `QLEN` (default `5000`)
- `TLEN` (default `5000`)

Example:

```bash
REPEAT=500 QLEN=8000 TLEN=8000 ./scripts/benchmark_compare.sh
```

## Test Example to Run

```bash
./tests/test_example.sh
```

What the example test checks:

- Scalar mode works on FASTA input.
- Auto mode returns same score as scalar.
- SIMD mode (if available) matches scalar score.
- The script prints traceback alignment output using `--show-alignment`.
