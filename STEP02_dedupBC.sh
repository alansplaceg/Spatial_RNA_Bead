#!/usr/bin/env bash
set -euo pipefail

for f in UDP8_extracted_barcodes.tsv; do
    out="${f%.tsv}_counted.tsv"

    echo "Processing $f"

    {
        printf "barcode_gene\tposID\tumi18\tn_reads\n"

        tail -n +2 "$f" \
        | LC_ALL=C sort --parallel=8 -S 50% \
        | awk -F'\t' '
            BEGIN {OFS="\t"}
            NR==1 {
                prev=$0
                n=1
                next
            }
            $0==prev {
                n++
                next
            }
            {
                print prev, n
                prev=$0
                n=1
            }
            END {
                if (NR) print prev, n
            }
        '
    } > "$out"

    echo "Wrote $out"
done
