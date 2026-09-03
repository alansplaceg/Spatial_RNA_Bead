import gzip, os
from concurrent.futures import ProcessPoolExecutor

FILES = [f"UDP{i}" for i in range(1, 10)]
MAX_WORKERS = 3

def reverse_complement(seq):
    complement = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(complement)[::-1]


def process_sample(sample):

    r1_file = f"{sample}_R1.fastq.gz"
    r2_file = f"{sample}_R2.fastq.gz"
    outfile = f"{sample}_extracted_barcodes.tsv"

    total = 0
    extracted = 0

    print(f"Processing {sample}", flush=True)

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


def main():

    files = [
        fq for fq in FILES
        if os.path.exists(fq)
    ]

    for fq in FILES:
        if not os.path.exists(fq):
            print("Missing:", fq)

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
