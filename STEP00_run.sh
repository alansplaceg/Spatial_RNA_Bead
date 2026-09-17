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
r1_only_count=0
for r1_file in "${r1_files[@]}"; do
    sample="${r1_file%_R1.fastq.gz}"
    r2_file="${sample}_R2.fastq.gz"

    if [[ -f "$r2_file" ]]; then
        ((pair_count += 1))
    else
        echo "No R2 found for $sample, will use R1-only anchor mode"
        ((r1_only_count += 1))
    fi
done

printf 'FASTQ directory: %s\n' "$TARGET_DIR"
printf 'Complete R1/R2 pairs: %d\n' "$pair_count"
printf 'R1-only samples: %d\n' "$r1_only_count"
printf 'Running STEP01...\n'

python3 "$SCRIPT_DIR/STEP01_getBCs.py"

for r1_file in "${r1_files[@]}"; do
    sample="${r1_file%_R1.fastq.gz}"

    printf '\nRunning downstream steps for %s...\n' "$sample"
    bash "$SCRIPT_DIR/STEP02_dedupBC.sh" "$sample"
    bash "$SCRIPT_DIR/STEP03_FindScaffold.sh" "$sample"
    bash "$SCRIPT_DIR/STEP04_getSameScaffold.sh" "$sample"
    bash "$SCRIPT_DIR/STEP05_annotateBCgenes.sh" "$sample"
done
