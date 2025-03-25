"""
Do not use this script. It does not handle multi-map reads. Use `samtools fastq` instead.
"""

import argparse
import gzip

import pysam


def bam_to_fastq(bam_file, output_fastq):
    """
    Convert a single-end BAM file to a gzipped FASTQ file using pysam.

    Args:
        bam_file (str): Path to the input BAM file.
        output_fastq (str): Path for the output gzipped FASTQ file.
    """
    with pysam.AlignmentFile(bam_file, "rb") as bam, gzip.open(output_fastq, "wt") as fq:
        for read in bam:
            fastq_entry = (
                f"@{read.query_name}\n{read.query_sequence}\n+\n{''.join(chr(q + 33) for q in read.query_qualities)}\n"
            )
            fq.write(fastq_entry)

    print(f"Conversion completed: {output_fastq}")


def main():
    parser = argparse.ArgumentParser(description="Convert single-end BAM to gzipped FASTQ using pysam.")
    parser.add_argument("bam_file", help="Input single-end BAM file")
    parser.add_argument("output_fastq", help="Output gzipped FASTQ file")
    args = parser.parse_args()

    bam_to_fastq(args.bam_file, args.output_fastq)


if __name__ == "__main__":
    main()
