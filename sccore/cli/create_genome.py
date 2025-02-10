#!usr/bin/env python3
"""
subset given chromsome from fasta and gtf
"""

import argparse
import pysam
import os


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gtf", required=True, help="gtf file")
    parser.add_argument("-f", "--fasta", required=True, help="fasta file")
    parser.add_argument(
        "-c", "--chrom", required=True, help="chromosome to extract. multple chrom are seperated by comma"
    )
    parser.add_argument("-o", "--species", required=True, help="use for output name")
    args = parser.parse_args()

    chrom = args.chrom.split(",")
    prefix = ".".join([args.species] + chrom)
    chrom = set(chrom)
    fa = pysam.FastxFile(args.fasta)
    outdir = prefix
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    out_fa = os.path.join(outdir, prefix + ".fasta")
    with open(out_fa, "w") as out:
        for entry in fa:
            name, seq = entry.name, entry.sequence
            if name in chrom:
                out.write(">{}\n{}\n".format(name, seq))

    out_gtf = os.path.join(outdir, prefix + ".gtf")
    gtf = open(args.gtf, "r")
    with open(out_gtf, "w") as out:
        for line in gtf:
            if line.split("\t")[0] in chrom:
                out.write(line)


if __name__ == "__main__":
    main()
