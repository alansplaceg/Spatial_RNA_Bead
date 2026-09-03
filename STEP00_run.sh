#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TARGET_DIR="$PWD"

shopt -s nullglob
r1_files=( *_R1.fastq.gz )

if (( ${#r1_files[@]} == 0 )); then
    echo "ERROR: No *_R1.fastq.gz files found in $TARGET_DIR" >&2
    exit 1
fi

pair_count=0
for r1_file in "${r1_files[@]}"; do
    sample="${r1_file%_R1.fastq.gz}"
    r2_file="${sample}_R2.fastq.gz"

    if [[ ! -f "$r2_file" ]]; then
        echo "ERROR: Missing pair for $r1_file: $r2_file" >&2
        exit 1
    fi

    ((pair_count += 1))
done

printf 'FASTQ directory: %s\n' "$TARGET_DIR"
printf 'Complete R1/R2 pairs: %d\n' "$pair_count"
printf 'Running STEP01...\n'

python3 "$SCRIPT_DIR/STEP01_getBCs.py"
