#!/usr/bin/env bash
set -euo pipefail

export LC_ALL=C

DIR="."
SAMPLE="${1:?Usage: STEP04_getSameScaffold.sh SAMPLE}"

IN="$DIR/${SAMPLE}_scaffold_matches.tsv"

RAW="$DIR/${SAMPLE}_scaffold_groups_raw.tsv"
FILTERED="$DIR/${SAMPLE}_scaffold_groups_filtered.tsv"

SORT_TMP="$DIR/sort_tmp"
KEEP="$SORT_TMP/${SAMPLE}_scaffolds_keep.txt"

THREADS=16
SORT_MEM="50%"

mkdir -p "$SORT_TMP"

[[ -f "$IN" ]] || {
    echo "ERROR: Missing $IN" >&2
    exit 1
}

echo "Input    : $IN"
echo "Raw      : $RAW"
echo "Filtered : $FILTERED"
echo


# ============================================================
# STEP 1
# RAW GROUPS
#
# Include UNIQUE + AMBIGUOUS.
# Exclude NO_MATCH.
#
# Collapse:
# scaffold_id + posID + barcode_gene
#
# Sum n_reads.
# ============================================================

echo "Building raw scaffold groups..."

{
    printf "scaffold_id\tposID\tbarcode_gene\tn_reads\n"

    awk -F'\t' -v OFS='\t' '
        NR > 1 && $7 != "NA" && ($8 == "UNIQUE" || $8 == "AMBIGUOUS") {
            print $7, $3, $2, $6
        }
    ' "$IN" \
    | sort \
        --parallel="$THREADS" \
        -S "$SORT_MEM" \
        -T "$SORT_TMP" \
        -t $'\t' \
        -k1,1n \
        -k2,2 \
        -k3,3 \
    | awk -F'\t' -v OFS='\t' '
        NR == 1 {
            scaffold=$1
            posid=$2
            barcode=$3
            reads=$4
            next
        }

        $1 == scaffold && $2 == posid && $3 == barcode {
            reads += $4
            next
        }

        {
            print scaffold, posid, barcode, reads

            scaffold=$1
            posid=$2
            barcode=$3
            reads=$4
        }

        END {
            if (NR > 0)
                print scaffold, posid, barcode, reads
        }
    '

} > "$RAW"

echo "Wrote: $RAW"


# ============================================================
# STEP 2
# FIND QUALIFYING SCAFFOLDS
#
# Require:
#
# >= 2 grouped rows
# >= 2 distinct posIDs
# >= 2 distinct barcode genes
# ============================================================

echo
echo "Finding qualifying scaffolds..."

awk -F'\t' '
    NR == 1 { next }

    {
        s=$1

        rows[s]++
        pos[s SUBSEP $2]=1
        bc[s SUBSEP $3]=1
    }

    END {
        for (x in pos) {
            split(x,a,SUBSEP)
            npos[a[1]]++
        }

        for (x in bc) {
            split(x,a,SUBSEP)
            nbc[a[1]]++
        }

        for (s in rows) {
            if (rows[s] >= 2 && npos[s] >= 2 && nbc[s] >= 2)
                print s
        }
    }
' "$RAW" \
| sort -n \
> "$KEEP"


# ============================================================
# STEP 3
# FILTER RAW GROUPS
#
# Keep ALL rows from qualifying scaffold IDs.
# ============================================================

awk -F'\t' -v OFS='\t' '
    NR == FNR {
        keep[$1]=1
        next
    }

    FNR == 1 {
        print
        next
    }

    ($1 in keep) {
        print
    }
' "$KEEP" "$RAW" \
> "$FILTERED"

rm -f "$KEEP"


# ============================================================
# SUMMARY
# ============================================================

echo
echo "RAW groups:"

awk -F'\t' '
    NR > 1 {
        scaffolds[$1]=1
        rows++
        reads += $4
    }

    END {
        ns=0
        for (x in scaffolds)
            ns++

        print "scaffolds=" ns, "rows=" rows, "reads=" reads
    }
' "$RAW"


echo
echo "FILTERED groups:"

awk -F'\t' '
    NR > 1 {
        scaffolds[$1]=1
        rows++
        reads += $4
    }

    END {
        ns=0
        for (x in scaffolds)
            ns++

        print "scaffolds=" ns, "rows=" rows, "reads=" reads
    }
' "$FILTERED"


echo
echo "Outputs:"
echo "  Raw:      $RAW"
echo "  Filtered: $FILTERED"
