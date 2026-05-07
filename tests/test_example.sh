#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

extract_value() {
    local line="$1"
    local key="$2"
    awk -v key="$key" '{
        for (i = 1; i <= NF; ++i) {
            split($i, kv, "=")
            if (kv[1] == key) {
                print kv[2]
                exit
            }
        }
    }' <<< "$line"
}

make >/dev/null

scalar_fasta="$(./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode scalar)"
scalar_score="$(extract_value "$scalar_fasta" best_score)"

auto_line="$(./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode auto)"
auto_score="$(extract_value "$auto_line" best_score)"

if [[ "$auto_score" != "$scalar_score" ]]; then
    echo "auto mode mismatch with scalar on FASTA input"
    exit 1
fi

if simd_line="$(./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode simd 2>/dev/null)"; then
    simd_score="$(extract_value "$simd_line" best_score)"
    if [[ "$simd_score" != "$scalar_score" ]]; then
        echo "SIMD score ($simd_score) does not match scalar score ($scalar_score)"
        exit 1
    fi
else
    echo "SIMD unavailable on this machine; explicit SIMD check skipped"
fi

echo "Example test passed"
echo "$scalar_fasta"
echo "Alignment output:"
./build/water --query tests/data/query.fasta --target tests/data/target.fasta --mode auto --show-alignment
