import gzip
import glob
import os
from concurrent.futures import ProcessPoolExecutor

MAX_WORKERS = 3

LINKER = "TTAGTGAGTTTGAGTTTGTT"
HYB    = "CGTTTCTGCTCTGTTCCCAG"

LINKER_MM = 2
HYB_MM = 2


def reverse_complement(seq):
    complement = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(complement)[::-1]


def mm_ok(a, b, k):
    mm = 0

    for x, y in zip(a, b):
        if x != y:
            mm += 1
            if mm > k:
                return False

    return True


def find_anchor(seq):
    """
    Structure (R1-only mode):

    barcode_gene | Linker | Hyb | posID | UMI18
         8           20      20      6      18

    Allows:
      <=2 mismatches in Linker
      <=2 mismatches in Hyb
    """

    for p in range(8, len(seq) - 64 + 1):

        if not mm_ok(
            seq[p:p+20],
            LINKER,
            LINKER_MM
        ):
            continue

        if mm_ok(
            seq[p+20:p+40],
            HYB,
            HYB_MM
        ):
            return p

    return None


def process_sample_r1_only(sample, r1_file):

    outfile = f"{sample}_extracted_barcodes.tsv"

    total = 0
    extracted = 0

    with gzip.open(r1_file, "rt") as f, open(outfile, "w") as out:

        out.write(
            "barcode_gene\tposID\tumi18\n"
        )

        while True:

            if not f.readline():
                break

            seq = f.readline().strip().upper()

            f.readline()

            if not f.readline():
                break

            total += 1

            p = find_anchor(seq)

            if p is None:
                continue

            barcode_gene = seq[p-8:p]
            posID = seq[p+40:p+46]
            umi18 = seq[p+46:p+64]

            if (
                len(barcode_gene) != 8
                or len(posID) != 6
                or len(umi18) != 18
            ):
                continue

            # Stream immediately to this sample's TSV
            out.write(
                f"{barcode_gene}\t"
                f"{posID}\t"
                f"{umi18}\n"
            )

            extracted += 1

    return sample, outfile, total, extracted


def process_sample_r1_r2(sample, r1_file, r2_file):

    outfile = f"{sample}_extracted_barcodes.tsv"

    total = 0
    extracted = 0

    with (
        gzip.open(r1_file, "rt") as r1,
        gzip.open(r2_file, "rt") as r2,
        open(outfile, "w") as out,
    ):

        out.write(
            "barcode_gene\tposID\tumi18\n"
        )

        while True:

            r1_header = r1.readline()
            r2_header = r2.readline()

            if not r1_header or not r2_header:
                break

            r1_seq = r1.readline().strip().upper()
            r2_seq = r2.readline().strip().upper()

            r1.readline()
            r2.readline()

            r1_quality = r1.readline()
            r2_quality = r2.readline()

            if not r1_quality or not r2_quality:
                break

            total += 1

            if len(r1_seq) < 24 or len(r2_seq) < 8:
                continue

            posID = r1_seq[:6]
            umi18 = r1_seq[6:24]
            barcode_gene = reverse_complement(r2_seq[:8])

            # Stream immediately to this sample's TSV
            out.write(
                f"{barcode_gene}\t"
                f"{posID}\t"
                f"{umi18}\n"
            )

            extracted += 1

    return sample, outfile, total, extracted


def process_sample(files):

    sample, r1_file, r2_file = files

    print(f"Processing {sample}", flush=True)

    if r2_file is None:
        return process_sample_r1_only(sample, r1_file)

    return process_sample_r1_r2(sample, r1_file, r2_file)


def main():

    files = []

    for r1_file in sorted(glob.glob("*_R1.fastq.gz")):
        sample = r1_file.removesuffix("_R1.fastq.gz")
        r2_file = f"{sample}_R2.fastq.gz"

        if os.path.exists(r2_file):
            files.append((sample, r1_file, r2_file))
        else:
            print(f"No R2 found for {sample}, using R1-only anchor mode")
            files.append((sample, r1_file, None))

    with ProcessPoolExecutor(
        max_workers=MAX_WORKERS
    ) as pool:

        results = list(
            pool.map(process_sample, files)
        )

    print()

    for sample, outfile, total, extracted in results:

        pct = (
            100 * extracted / total
            if total else 0
        )

        print(
            f"{sample}: "
            f"{extracted:,}/{total:,} "
            f"({pct:.2f}%) -> {outfile}"
        )


if __name__ == "__main__":
    main()
