#!/usr/bin/env bash
set -euo pipefail

sample="${1:?Usage: STEP05_annotateBCgenes.sh SAMPLE}"
IN="${sample}_scaffold_groups_filtered.tsv"
OUT="${sample}_scaffold_groups_filtered_annotated.tsv"

[[ -f "$IN" ]] || {
    echo "ERROR: Missing $IN" >&2
    exit 1
}

awk -F'\t' -v OFS='\t' '

function hd(a,b,    i,n) {
    n=0

    if (length(a) != length(b))
        return 99

    for (i=1; i<=length(a); i++)
        if (substr(a,i,1) != substr(b,i,1))
            n++

    return n
}

BEGIN {

    # posID annotation
    posname["TAGTTG"]="PosID1"
    posname["TACAAG"]="PosID2"
    posname["AGGAAT"]="PosID3"
    posname["CTTTTG"]="PosID4"


    # 8-nt barcode whitelist
    probe["AAGGTCCA"]="GAPDH Probe 1"
    probe["AGAGTACC"]="GAPDH Probe 2"
    probe["CTGTACCA"]="GAPDH Probe 3"
    probe["GCAGACAA"]="GAPDH Probe 4"
    probe["CCAGTCTT"]="GAPDH Probe 5"
    probe["TACCAGGT"]="GAPDH Probe 6"
    probe["TCTCGATC"]="GAPDH Probe 7"
    probe["GTCAAGTG"]="GAPDH Probe 8"
    probe["GCTAAGGT"]="GAPDH Probe 9"
    probe["TGCTCCAA"]="GAPDH Probe 10"
    probe["TCGATTGG"]="GAPDH Probe 11"
    probe["TACCACCT"]="GAPDH Probe 12"
    probe["AGCTGGAT"]="GAPDH Probe 13"
    probe["TACACTCC"]="GAPDH Probe 14"
    probe["TTGACGGA"]="GAPDH Probe 15"
    probe["AGGAACTC"]="GAPDH Probe 16"
    probe["CCTTGTCA"]="GAPDH Probe 17"
    probe["GTACGACT"]="GAPDH Probe 18"
    probe["GAGGACTT"]="GAPDH Probe 19"
    probe["GCTTCTTG"]="GAPDH Probe 20"
    probe["ACTCCTCT"]="GAPDH Probe 21"
    probe["ATGCCTAC"]="GAPDH Probe 22"
    probe["GATGTGTG"]="GAPDH Probe 23"
    probe["TGGATGCT"]="GAPDH Probe 24"

    probe["ACCTGAGT"]="ACTB Probe 1"
    probe["GTGACTAC"]="B2M Probe 1"
    probe["TCAAGTCC"]="PPIA Probe 1"
    probe["CGTTAAGC"]="UBC Probe 1"
    probe["ATCCGTGA"]="EEF2 Probe 1"


    # Full Variable Region B
    var["AAGGTCCA"]="AAGGTCCAT"
    var["AGAGTACC"]="AGAGTACCT"
    var["CTGTACCA"]="CTGTACCAT"
    var["GCAGACAA"]="GCAGACAAT"
    var["CCAGTCTT"]="CCAGTCTTT"
    var["TACCAGGT"]="TTACCAGGTT"
    var["TCTCGATC"]="TCTCGATCT"
    var["GTCAAGTG"]="GTCAAGTGT"
    var["GCTAAGGT"]="GCTAAGGTT"
    var["TGCTCCAA"]="TGCTCCAAT"
    var["TCGATTGG"]="TCGATTGGT"
    var["TACCACCT"]="TACCACCTT"
    var["AGCTGGAT"]="AGCTGGATT"
    var["TACACTCC"]="TACACTCCT"
    var["TTGACGGA"]="TTGACGGAT"
    var["AGGAACTC"]="AGGAACTCT"
    var["CCTTGTCA"]="CCTTGTCAT"
    var["GTACGACT"]="GTACGACTT"
    var["GAGGACTT"]="GAGGACTTT"
    var["GCTTCTTG"]="GCTTCTTGT"
    var["ACTCCTCT"]="ACTCCTCTT"
    var["ATGCCTAC"]="ATGCCTACT"
    var["GATGTGTG"]="GATGTGTGT"
    var["TGGATGCT"]="TGGATGCTT"

    var["ACCTGAGT"]="ACCTGAGTT"
    var["GTGACTAC"]="GTGACTACT"
    var["TCAAGTCC"]="TCAAGTCCT"
    var["CGTTAAGC"]="CGTTAAGCT"
    var["ATCCGTGA"]="ATCCGTGAT"
}

NR == 1 {
    print \
        "scaffold_id", \
        "posID", \
        "posID_name", \
        "barcode_gene", \
        "corrected_barcode", \
        "probe_name", \
        "variable_region_B", \
        "barcode_mismatches", \
        "n_reads"

    next
}

{
    pos=$2
    bc=$3

    # Require valid posID
    if (!(pos in posname))
        next

    pname=posname[pos]


    # Exact barcode
    if (bc in probe) {
        best=bc
        mismatches=0
    }

    else {
        # Unique <=1 mismatch barcode correction
        best=""
        nhits=0

        for (x in probe) {
            if (hd(bc,x) <= 1) {
                nhits++
                best=x
            }
        }

        # Remove unannotated AND ambiguous 1-MM rows
        if (nhits != 1)
            next

        mismatches=1
    }


    print \
        $1, \
        pos, \
        pname, \
        bc, \
        best, \
        probe[best], \
        var[best], \
        mismatches, \
        $4
}

' "$IN" > "$OUT"


echo "Wrote: $OUT"

echo
echo "Annotation summary:"

awk -F'\t' '
NR > 1 {
    rows[$6]++
    reads[$6]+=$9
}
END {
    for (x in rows)
        print x, "rows=" rows[x], "reads=" reads[x]
}
' "$OUT" | sort
