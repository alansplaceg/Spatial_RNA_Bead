#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

DIR="."
SAMPLE="${1:?Usage: STEP03_FindScaffold.sh SAMPLE}"

SCAFFOLDS="/mnt/spatialdata/FastqData/1HBpreseq/Bub_barcodes_combined.tsv"
READS="$DIR/${SAMPLE}_extracted_barcodes_counted.tsv"

OUT="$DIR/${SAMPLE}_scaffold_matches.tsv"
INDEX="$DIR/Bub_position_index.tsv"
INDEX_META="$DIR/Bub_position_index.source"

SORT_TMP="$DIR/sort_tmp"
READ_KEYS="$SORT_TMP/${SAMPLE}_read_keys.tsv"

THREADS=16
SORT_MEM="50%"

mkdir -p "$SORT_TMP"

echo "Scaffolds : $SCAFFOLDS"
echo "Reads     : $READS"
echo "Output    : $OUT"
echo


# ============================================================
# CHECK INPUT FILES
# ============================================================

[[ -f "$SCAFFOLDS" ]] || {
    echo "ERROR: Missing $SCAFFOLDS" >&2
    exit 1
}

[[ -f "$READS" ]] || {
    echo "ERROR: Missing $READS" >&2
    exit 1
}


# ============================================================
# STEP 1
# BUILD REUSABLE SCAFFOLD INDEX
#
# Bub1 -> position 1
# Bub2 -> position 2
# Bub3 -> position 3
# Bub4 -> position 4
#
# Index format:
#
# key                scaffold_id(s)     status
#
# 1:ACGT...          1234               UNIQUE
# 2:TGCA...          100,500,900        AMBIGUOUS
#
# For ambiguous barcodes, retain ALL matching scaffold IDs.
# ============================================================

SCAFFOLD_SIGNATURE="$(stat -c '%s:%Y' "$SCAFFOLDS")"
CACHED_SIGNATURE=""
if [[ -s "$INDEX_META" ]]; then
    CACHED_SIGNATURE="$(cat "$INDEX_META")"
fi

if [[ ! -s "$INDEX" || "$CACHED_SIGNATURE" != "$SCAFFOLDS"$'\n'"$SCAFFOLD_SIGNATURE" ]]; then

    echo "Building scaffold index..."

    awk -F'\t' -v OFS='\t' '
        NR == 1 {
            next
        }

        {
            scaffold_id = NR - 1

            print "1:" $1, scaffold_id
            print "2:" $2, scaffold_id
            print "3:" $3, scaffold_id
            print "4:" $4, scaffold_id
        }
    ' "$SCAFFOLDS" \
    | sort \
        --parallel="$THREADS" \
        -S "$SORT_MEM" \
        -T "$SORT_TMP" \
        -t $'\t' \
        -k1,1 \
        -k2,2n \
    | awk -F'\t' -v OFS='\t' '

        function flush_key() {

            if (key == "")
                return

            if (n == 1)
                print key, ids, "UNIQUE"
            else
                print key, ids, "AMBIGUOUS"
        }

        {
            if ($1 != key) {

                flush_key()

                key = $1
                ids = $2
                n = 1
            }
            else {

                ids = ids "," $2
                n++
            }
        }

        END {
            flush_key()
        }

    ' > "${INDEX}.tmp"

    mv "${INDEX}.tmp" "$INDEX"
    printf '%s\n%s\n' "$SCAFFOLDS" "$SCAFFOLD_SIGNATURE" > "${INDEX_META}.tmp"
    mv "${INDEX_META}.tmp" "$INDEX_META"

    echo "Index created:"
    echo "$INDEX"

else

    echo "Using existing scaffold index:"
    echo "$INDEX"

fi


# ============================================================
# STEP 2
# PREPARE READ LOOKUP KEYS
#
# Input columns:
#
# barcode_gene  posID  umi18  n_reads
#
# Reverse-complement umi18 before scaffold lookup.
#
# posID mapping in READ orientation:
#
# TAGTTG -> PosID1 -> Bub1
# TACAAG -> PosID2 -> Bub2
# AGGAAT -> PosID3 -> Bub3
# CTTTTG -> PosID4 -> Bub4
#
# Lookup key:
#
# position:reverse_complement(umi18)
# ============================================================

echo
echo "Preparing read keys..."

awk -F'\t' -v OFS='\t' '

    function revcomp(seq,    i,b,out) {

        out = ""

        for (i = length(seq); i >= 1; i--) {

            b = substr(seq, i, 1)

            if      (b == "A") b = "T"
            else if (b == "T") b = "A"
            else if (b == "C") b = "G"
            else if (b == "G") b = "C"
            else if (b == "N") b = "N"

            out = out b
        }

        return out
    }


    NR == 1 {
        next
    }


    {
        row_id = NR - 1

        barcode_gene = $1
        posID = $2
        umi18 = $3
        n_reads = $4

        umi18_rc = revcomp(umi18)


        if (posID == "TAGTTG") {
            key = "1:" umi18_rc
        }
        else if (posID == "TACAAG") {
            key = "2:" umi18_rc
        }
        else if (posID == "AGGAAT") {
            key = "3:" umi18_rc
        }
        else if (posID == "CTTTTG") {
            key = "4:" umi18_rc
        }
        else {
            key = "X:" umi18_rc
        }


        # key
        # row_id
        # barcode_gene
        # posID
        # umi18
        # umi18_rc
        # n_reads

        print \
            key, \
            row_id, \
            barcode_gene, \
            posID, \
            umi18, \
            umi18_rc, \
            n_reads
    }

' "$READS" \
| sort \
    --parallel="$THREADS" \
    -S "$SORT_MEM" \
    -T "$SORT_TMP" \
    -t $'\t' \
    -k1,1 \
> "$READ_KEYS"


# ============================================================
# STEP 3
# MATCH READS TO SCAFFOLDS
#
# UNIQUE
#   one matching scaffold
#   -> one output row
#
# AMBIGUOUS
#   multiple matching scaffolds
#   -> one output row PER matching scaffold
#
# NO_MATCH
#   no matching scaffold
#   -> scaffold_id = NA
#
# Output columns:
#
# row_id
# barcode_gene
# posID
# umi18
# umi18_rc
# n_reads
# scaffold_id
# status
# ============================================================

echo "Matching reads to scaffolds..."

{
    printf "row_id\tbarcode_gene\tposID\tumi18\tumi18_rc\tn_reads\tscaffold_id\tstatus\n"

    join \
        -t $'\t' \
        -a 2 \
        -e "NA" \
        -o '2.2,2.3,2.4,2.5,2.6,2.7,1.2,1.3' \
        "$INDEX" \
        "$READ_KEYS" \
    | awk -F'\t' -v OFS='\t' '

        {
            row_id       = $1
            barcode_gene = $2
            posID        = $3
            umi18        = $4
            umi18_rc     = $5
            n_reads      = $6
            scaffold_ids = $7
            status       = $8


            # --------------------------------------------
            # NO MATCH
            # --------------------------------------------

            if (scaffold_ids == "NA") {

                print \
                    row_id, \
                    barcode_gene, \
                    posID, \
                    umi18, \
                    umi18_rc, \
                    n_reads, \
                    "NA", \
                    "NO_MATCH"

                next
            }


            # --------------------------------------------
            # UNIQUE MATCH
            # --------------------------------------------

            if (status == "UNIQUE") {

                print \
                    row_id, \
                    barcode_gene, \
                    posID, \
                    umi18, \
                    umi18_rc, \
                    n_reads, \
                    scaffold_ids, \
                    "UNIQUE"

                next
            }


            # --------------------------------------------
            # AMBIGUOUS MATCH
            #
            # Example:
            #
            # scaffold_ids = 123,456,789
            #
            # becomes:
            #
            # ... 123 AMBIGUOUS
            # ... 456 AMBIGUOUS
            # ... 789 AMBIGUOUS
            # --------------------------------------------

            if (status == "AMBIGUOUS") {

                n_ids = split(scaffold_ids, ids, ",")

                for (i = 1; i <= n_ids; i++) {

                    print \
                        row_id, \
                        barcode_gene, \
                        posID, \
                        umi18, \
                        umi18_rc, \
                        n_reads, \
                        ids[i], \
                        "AMBIGUOUS"
                }

                next
            }


            # --------------------------------------------
            # FALLBACK
            # --------------------------------------------

            print \
                row_id, \
                barcode_gene, \
                posID, \
                umi18, \
                umi18_rc, \
                n_reads, \
                "NA", \
                "NO_MATCH"
        }

    '

} > "$OUT"


# ============================================================
# REMOVE TEMP READ KEY FILE
# ============================================================

rm -f "$READ_KEYS"


# ============================================================
# SUMMARY
# ============================================================

echo
echo "Finished."
echo "Output:"
echo "$OUT"
echo


# ============================================================
# SUMMARY BY STATUS
#
# NOTE:
# AMBIGUOUS reads are expanded to multiple scaffold rows.
# Therefore AMBIGUOUS read totals here may exceed the number
# of original sequencing reads.
# ============================================================

echo "Summary by status:"

awk -F'\t' '
    NR > 1 {
        rows[$8]++
        reads[$8] += $6
    }

    END {
        for (x in rows)
            print x, "rows=" rows[x], "reads=" reads[x]
    }
' "$OUT" \
| sort


# ============================================================
# SUMMARY BY posID
# ============================================================

echo
echo "Summary by posID:"

awk -F'\t' '
    NR > 1 {
        rows[$3]++
        reads[$3] += $6
    }

    END {
        for (x in rows)
            print x, "rows=" rows[x], "reads=" reads[x]
    }
' "$OUT" \
| sort


# ============================================================
# UNIQUE + AMBIGUOUS MATCHES
# ============================================================

echo
echo "Matched scaffold rows (UNIQUE + AMBIGUOUS):"

awk -F'\t' '
    NR > 1 && ($8 == "UNIQUE" || $8 == "AMBIGUOUS") {
        rows++
        reads += $6
    }

    END {
        print "rows=" rows, "reads=" reads
    }
' "$OUT"


# ============================================================
# UNIQUE ONLY
# ============================================================

echo
echo "Unique matches:"

awk -F'\t' '
    NR > 1 && $8 == "UNIQUE" {
        rows++
        reads += $6
    }

    END {
        print "rows=" rows, "reads=" reads
    }
' "$OUT"


# ============================================================
# AMBIGUOUS ONLY
# ============================================================

echo
echo "Ambiguous scaffold rows:"

awk -F'\t' '
    NR > 1 && $8 == "AMBIGUOUS" {
        rows++
        reads += $6
    }

    END {
        print "rows=" rows, "reads=" reads
    }
' "$OUT"
