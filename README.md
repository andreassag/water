> [!NOTE]
> **Archival Notice**: This repository is complete, archived, and in read-only mode. It is preserved for educational and reference purposes; no further active maintenance, feature additions, or pull requests will be processed.

# SIMD Smith-Waterman in C

This project implements local sequence alignment (Smith-Waterman) with affine gap penalties in two versions:

- Scalar dynamic programming baseline.
- SIMD-optimized striped implementation with runtime-dispatched AVX2 or SSE2 intrinsics.

Both versions support the same scoring parameters and can be selected from the CLI.

## Features

- Affine gap penalty model:
  - Gap-open penalty
  - Gap-extend penalty
- Input format support: FASTA and FASTQ files.
- CLI alignment mode selection:
  - `scalar`
  - `simd`
  - `auto` (detects AVX2 or SSE2 at runtime; falls back to scalar if unavailable)
- Benchmark script to compare scalar vs SIMD throughput.
- Runtime SIMD backend reporting (`simd_backend`, `active_backend`) in CLI output.
- Optional full visual alignment traceback (`--show-alignment`).

## Prerequisites & Build Instructions

### Prerequisites

- A C11-compliant C compiler (`gcc` or `clang`).
- GNU `make`.
- Linux x86_64 system (AVX2 or SSE2 instruction set support recommended for SIMD acceleration).

### Building

Build the optimized executable:

```bash
make
```

Binary output:
- `build/water`

You can customize the compiler or flags by passing variables to `make`:

```bash
make CC=gcc CFLAGS="-O3 -march=native"
# Or with clang (if installed):
# make CC=clang CFLAGS="-O3 -march=native"
```

### Running Tests

Run the integrated test suite:

```bash
make test
```

This compiles `build/water` if necessary and executes [tests/test_example.sh](file:///home/andreassag/repo/water/tests/test_example.sh), validating FASTA parsing, score parity between scalar and SIMD modes, and traceback visualization.

### Cleaning Build Artifacts

To remove compiled object files and binaries:

```bash
make clean
```

## CLI Usage

### Command Syntax

```bash
./build/water --query <query.fa|fq> --target <target.fa|fq> [options]
```

### Options

| Option | Description | Default |
| :--- | :--- | :--- |
| `-q, --query <path>` | Path to query FASTA or FASTQ file | *(Required)* |
| `-t, --target <path>` | Path to target FASTA or FASTQ file | *(Required)* |
| `-m, --mode <mode>` | Alignment mode: `auto`, `scalar`, or `simd` | `auto` |
| `--match <int>` | Score reward for matching bases | `2` |
| `--mismatch <int>` | Penalty score for mismatching bases | `-2` |
| `--gap-open <int>` | Penalty score for opening a gap (`> 0`) | `3` |
| `--gap-extend <int>`| Penalty score for extending a gap (`> 0`) | `1` |
| `-r, --repeat <N>` | Repeat alignment computation `N` times for benchmarking | `1` |
| `-a, --all` | Align all pairwise combinations across multi-record files | Off |
| `-A, --show-alignment` | Print visual traceback alignment for the first query-target pair | Off |
| `-h, --help` | Display command usage and option descriptions | |

### Examples

**1. Basic local alignment (auto dispatch):**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode auto
```

**2. Explicit scalar mode (portable baseline):**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode scalar
```

**3. Explicit SIMD mode (vectorized acceleration):**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode simd
```

**4. Custom scoring parameters (e.g., match +3, mismatch -3, gap open 4, gap extend 2):**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --match 3 --mismatch -3 --gap-open 4 --gap-extend 2
```

**5. Benchmarking execution loops (repeat 10 iterations):**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --repeat 10
```

**6. Pairwise alignment of multi-record FASTA/FASTQ files:**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --all
```

**7. Visualizing the alignment traceback:**

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --show-alignment
```

## Sample Alignment Output

Running `water` with `--show-alignment` produces machine-parsable key-value metrics followed by the reconstructed sequence alignment:

```bash
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --show-alignment
```

### Output Block

```text
requested_mode=auto selected_mode=simd simd_backend=avx2 active_backend=avx2 alignments=1 repeat=1 best_score=785 score_sum=785 avg_ms=1.342742 gcups=1.335983 query_count=1 target_count=1
alignment_query=NZ_LR134336.1:c1928227-1926908 alignment_target=NZ_CP077224.1:51427-52785 alignment_score=785 query_range=0:1310 target_range=0:1349
query_aln=ATGAGTTTATTAGGATTTGGAATAGTT--ATTCTT--CA--TATT-T-----TTATT--TTAACA-ATA-AGC-TATGAATTAATT-TGTAAAGTAACTAAAAAAGAAAGATAT-TTATTTATAGCTTGTA-TATTT-GTATATCAGTCT---GTAGCTGAATT-TTTCTTTGACT-TTATT-TCATT---AGCTGGGTTAGGGAT-TGA-AA--AG----TT----TATA-TTTCCT-TTTATTCTA-----TATACTTCTA-TAGTAGTT--TT-AAAGA----AGTATGATAAGCATAAAGGAA--TG--TTTATCAGTTTACTATTGTCTTTATTATATCATAGT--ACAC-----ATACTTTTTTATCAGTAACCCTATCCTCTATAACAGGTGATGCTTTTGTGATAGAACATGAAGTGATGTTTTA-TTTAGGAGTACTATTACTAACATATTTTATTATTAAAAA--------AATTATT--AC--TTA-TTTCCAT-TTAGAACTTACCTATTTTGATAAAGATTATCTGTATCCTTTTTTAAAAAAAGTAATAGTTGCTTTCTTTGGT--TTACATATCTTATTGTTTATTTCAGATATAG-TAAG-TACAACTCATCATTTAAATAGTTTCGGAAGTATTTTATCAACCATTG--TTT-TTATTTGCTTGT-TA-TTGATTTTCTTTGCAATGAATTCTC----ACAAGGTTCAAATTGAAAAAGAGATT--GCTCTAAAACAAAAGAAATTTGAACAAAAACATTTACAGACTTATACTGATGAAATTGTGGAGT-TATATAATGAAATTCGCGGTTTTCGTCATGACTATGCAGGGATGCTTGTTAGCATGCAAATGGCGATTGACAGTGGAGATTTACAAGAAATTAACAGGGT-TTACAATGAAGTCTT-AGTAAAAGCAAATCAGAAATT------GAGGTCAGAAAAATATACTTACTTTGATTTAAATAATATAGAAGATTCTGCTTTACGAAGT-T--TAATTGCTCAATCGATTGTCTATGCACGAAAAAATGATGTAGAGTTTACATTGGAAGTAAAGGA---TAT-----CATTA-CTA-G--ACTGTCAATAGATTTACTT-GA-T--CTTGTTCGTATCATGAGCATTTTATTAAACAATGCTGTTGAAGGTGCTGCAGATAGTTATTTGAAACAAATGGAAGTTGCAGTCATTA--AAA--TGGATTTTGAAACAGTTATAGTTATCCAGAATTCATG-CAAA---ATTACTATGACGCCTTCAGAGGATCTATTTGCCTTAGGTTTCTCTACCAAGGGAAGAAATAGAGGTCTAGGACTTAATAATGTCAAGGAGATTTTAGATAAGTATGACAACATTATTTTAGAAACAGAGATGGAAGACAACACA-TTTAGACAAATTATTAGATTTAAGAG---GGAA
match_aln=|||| ||| |.|.|||    ||||.||  |||.||  ||  |||| |     ||.||  |||..| ||| ..| ||.|||  |||| ||.|.|||     |||||||.|||||| ..||||||.|.|| || |.||| |||| |  ||.|   |.|..|.|||| ||||||..|.| ||||| ..|||   ||.|   ||.|.||| |.| ||  ||    ||    |||| |||||| |||||||.|     ||| .||||| .|.||.||  || |||||    ||||||||||  .||...|||  ||  |||||| |||   |.||.||  |||| |||| | ||  |.||     |.|||||||||||.||.|...|.||.||.||||||||.|||..||||||.|.|.|||||.|.|..||.||||| .|||..|.||.|.|| .||.|.|..|||||||||..|||        |.||.||  ||  ||| |||..|| |||.||...|||.||||  |||||.||.||||.|||..|.||    ||.||||.|||         |.|.||  || |||||||||||...||||||..|| |.| |.|| .|.|||.|| |||.|.||||||||.|..|||||..|..|.||.||.|  ||| ||||.|..|||| .| ||.|||||.|||  ||.|..|.|||    |||| ..|.|||..|.||||||.|||  ||    ||| |||||.||   |||||||.|||..|||||..||||||.|||||||||||.| || |.|||||.||||||||.||||||||||||||.|||||.||.|||.|..|.||..|.||||..|..|||.||...||.|||.|..|||||.|..|.||..| ||.| ||.|.||.|| ||| .||||.|||    |.||      ||  |||||..|.|||||||..|||||.|||||||||.|..||||..||||.||..||||| |  ||||   ||||.|.||..|..|.||.|||.|...||..||||||.|...|||.|||.|.|||||   |||     ||||| ||| |  |||.|   ||||||||.|| || |  ||.||  ||||  |  ||       |.||.|||||.||||||||.||.||.||.|||...|..|||..||||.|.||    |||.|||  |||  |.|||...||||..|||.|.||.||.||||||||..| ||||   .||| |||.|   .||.|||.||..|.|.||...||||.|||||.||.||.|||..||||.||||..||||..|.|||||||||||.||.|||.|.||.||.|||||..|..|||||.|.||||||||.||.||| ||||...| ||||.||||.||.|.||||||||.||   ||||
target_aln=ATGA-TTT-TGAAGAT----AATAATTGCATTTTTTCCAACTATTATACAAGTTGTTTCTTATGAGATATTTCATAAGAA--AATTATGCAGAGT-----AAAAAGATAGATATAAAATTTATTGGTT-TATTCTTTAGTAT-T--GTTTTTGGAATTTCAATTATTTCTTCAATTCTTATTGAAATTCCAAGTT---TTGGAGATGTAATAAGTAGTATTTTCCACTATATTTTCCTATTTATTCAACCCCTTAT-TTTCTACAAATATTTTCTTAAAAGAAAAGAGTATGATAA--TTACTTGAATTTGTTTTTATC-GTT---TGTTATC--TATT-TATC-T-GTGGAGACTTCCGAAACTTTTTTATCTGTCATTATTTCTTCAATAACAGGAGATTATTTTGTAAAACAACATTACGATATTTTTTACATTATAATTAATCTT-TTATCTTTGTTTATTATTTTAAAGGTTGTTGATTTTTTTGACTTTTATTTTGAATATTATAAAGAACCAATTT--ATAAAAATGATCTATATAATGTT----AATAAGTCATA---------TATAGTAATT-CATATCTTATTAAATATTTCTCAT-TGGTTTAGCGAAAACGCA-CATCTTAATAGTTTTGCTAGTATGATTGCCACTATAGGATTTATTATGTTTTTGTCGACTTTATTTTATTT--AAAGTCTGCTCGTGAACAA-CATGAAAAAGCAAAAGAAATTCAGC----AAA-AAAAGGAA---GAACAAAGACAACTACAGCTTTATACAGATGAAATTGTAG-GTCTTTATAACGAAATTCGAGGTTTTCGTCATGATTATGCTGGTATGTTAATAAGTTTACAAACAGGAATTAACTCAGGTGATATGAAAGAAGTAGAGAGAATCTTTC-ATAATGTTTTAAGT-CAAGCTAAT----ATTTCCCTTCGA--TCAGATGACTATACTTTTTTTGAGTTAAATAATGTTCAAGACACTGCCTTGAGAAGTGTGCTAAT---TCAAACTATATTTAAGGCTCGAGAGTGTGGAGTAGAGATCGTATTTGAAATGAAGGATGTTATTGAGACATTACCTATGAAACTTT---TAGATTTAGTTCGAGTAGCTAGT--GTAT--T--GC-------TCAATAATGCGGTTGAAGGGGCAGCTGAAAGTCCCTCAAAAACAATGAACGT----GTCTTTAGTAAAGCTAGATAAAGAAATTGTTTTTGTAATTCAGAATTCTCGTCAAAGCCGTTA-TATTA---ATTTAGAAGAAATTTATGAAGTAGGCTTCTCAACAAAAGGAGAAAATCGAGGATTAGGTTTAAATAATGTCAAAGAAATTATTGACAAATATGATGAAGTTATTCTGGAAACAGATATAGAA-ACAAATTATTTTATACAAGTTGTAAGATTTAAAAGAAAGGAA
```

### Field Explanations

#### Summary Metrics Line
- `requested_mode`: The mode passed via `--mode` (`auto`, `scalar`, or `simd`).
- `selected_mode`: The alignment engine chosen (`scalar` or `simd`).
- `simd_backend`: The highest SIMD capability detected on the host CPU (`avx2`, `sse2`, or `none`).
- `active_backend`: The specific computational kernel executed (`avx2`, `sse2`, or `scalar`).
- `alignments`: Number of sequence pairs evaluated.
- `repeat`: Benchmark iterations completed per pair.
- `best_score`: The highest local alignment score discovered.
- `score_sum`: The accumulated alignment score across all runs.
- `avg_ms`: The mean runtime in milliseconds per benchmark repetition.
- `gcups`: Performance metric in **Giga Cell Updates Per Second**:
  $$\text{GCUPS} = \frac{\text{Total DP Matrix Cells} \times \text{Repeat Count}}{\text{Elapsed Seconds} \times 10^9}$$
  *(For a single sequence pair, $\text{Total DP Cells} = \text{Query Length} \times \text{Target Length}$.)*
- `query_count` / `target_count`: Number of FASTA/FASTQ records ingested from query and target files.

#### Alignment Details & Visual Traceback
- `alignment_query`: The FASTA/FASTQ identifier of the query sequence.
- `alignment_target`: The FASTA/FASTQ identifier of the target sequence.
- `alignment_score`: Optimal Smith-Waterman local alignment score for this pair.
- `query_range`: 0-based coordinate range `[start:end]` in the query sequence covered by the local alignment.
- `target_range`: 0-based coordinate range `[start:end]` in the target sequence covered by the local alignment.
- `query_aln`: Aligned query sequence (including `-` for gap insertions).
- `match_aln`: Visual match representation:
  - `|`: Exact base match.
  - `.`: Base mismatch or substitution.
  - ` `: Gap character (insertion/deletion).
- `target_aln`: Aligned target sequence (including `-` for gap deletions).

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
