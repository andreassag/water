#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

REPEAT="${REPEAT:-200}"
QLEN="${QLEN:-5000}"
TLEN="${TLEN:-5000}"

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

generate_fasta() {
    local out="$1"
    local name="$2"
    local len="$3"
    {
        echo ">$name"
        awk -v len="$len" 'BEGIN {
            srand();
            bases = "ACGT";
            for (i = 0; i < len; ++i) {
                idx = int(rand() * 4) + 1;
                printf "%s", substr(bases, idx, 1);
            }
            printf "\n";
        }'
    } > "$out"
}

make >/dev/null

tmp_dir="$(mktemp -d)"
trap 'rm -rf "$tmp_dir"' EXIT

query_file="$tmp_dir/query.fasta"
target_file="$tmp_dir/target.fasta"

generate_fasta "$query_file" query "$QLEN"
generate_fasta "$target_file" target "$TLEN"

scalar_line="$(./build/water --query "$query_file" --target "$target_file" --mode scalar --repeat "$REPEAT")"

simd_line=""
if simd_line="$(./build/water --query "$query_file" --target "$target_file" --mode simd --repeat "$REPEAT" 2>/dev/null)"; then
    scalar_ms="$(extract_value "$scalar_line" avg_ms)"
    simd_ms="$(extract_value "$simd_line" avg_ms)"

    speedup="$(awk -v s="$scalar_ms" -v v="$simd_ms" 'BEGIN { if (v == 0) { print "inf" } else { printf "%.2fx", s / v } }')"

    echo "$scalar_line"
    echo "$simd_line"
    echo "speedup=$speedup"
else
    echo "$scalar_line"
    echo "simd_status=unavailable"
fi
