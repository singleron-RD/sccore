import argparse
import gzip
from collections import defaultdict

import pysam
from sccore import parse_chemistry, utils

CHEMISTRYS = ["GEXSCOPE-V2", "bulk_rna-bulk_vdj_match"]
BULK_LINKER = "GTGGTATCAACGCAGAGT"


@utils.add_log
def run(args):
    fq1_list = args.fq1.split(",")
    fq2_list = args.fq2.split(",")
    auto_runner = parse_chemistry.AutoRNA(fq1_list)
    L18_set = parse_chemistry.create_mismatch_seqs(BULK_LINKER, max_mismatch=2)
    out_fq1 = {chemistry: gzip.open(f"{chemistry}_R1.fq.gz", mode="wt", compresslevel=1) for chemistry in CHEMISTRYS}
    out_fq2 = {chemistry: gzip.open(f"{chemistry}_R2.fq.gz", mode="wt", compresslevel=1) for chemistry in CHEMISTRYS}

    total = 0
    count_dict = defaultdict(int)

    for fq1_file, fq2_file in zip(fq1_list, fq2_list):
        with pysam.FastqFile(fq1_file) as f1, pysam.FastqFile(fq2_file) as f2:
            for read1, read2 in zip(f1, f2):
                total += 1
                chemistry = None
                if read1.sequence[:18] in L18_set:
                    chemistry = "bulk_rna-bulk_vdj_match"
                elif auto_runner.is_chemistry(read1.sequence, "GEXSCOPE-V2"):
                    chemistry = "GEXSCOPE-V2"
                if chemistry:
                    count_dict[chemistry] += 1
                    out_fq1[chemistry].write(str(read1) + "\n")
                    out_fq2[chemistry].write(str(read2) + "\n")

    print(f"total reads: {total}")
    for chemistry in CHEMISTRYS:
        print(f"{chemistry} reads: {count_dict[chemistry]}")
        out_fq1[chemistry].close()
        out_fq2[chemistry].close()


def main():
    parser = argparse.ArgumentParser(description="Split fastq file by chemistry")
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
