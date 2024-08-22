#!/usr/bin/env python

import sys
import pysam


def main():
    bam = sys.argv[1]
    total = 0
    diff_fh = open("diff_gene.tsv", "w")
    fh = pysam.AlignmentFile(bam, "rb")
    for read in fh:
        total += 1
        if total % 100000 == 0:
            print(f"Processed {total} reads")
        if not read.has_tag("XT"):
            continue
        gene = read.get_tag("XT")
        starsolo_gene = read.get_tag("GX")
        if gene != starsolo_gene:
            diff_fh.write(str(read) + "\n")

    print(f"Total reads: {total}")


if __name__ == "__main__":
    main()
