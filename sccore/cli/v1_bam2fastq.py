"""
Convert celescope V1.* BAM to fastq file
"""

import pysam
import gzip
import os

LINKER1 = "ATCCACGTGCTTGAGA"
LINKER2 = "TCAGCATGCGGCTACG"


def fastq_str(name, seq, qual):
    """return fastq read string"""
    return f"@{name}\n{seq}\n+\n{qual}\n"


def seg2fastq(segment: pysam.AlignedSegment, cb_len: int) -> tuple[str, str]:
    """
    C9L16C9L16C9L1U12
    """
    query_name = segment.query_name
    attr = query_name.split("_")

    cb, umi = attr[0], attr[1]
    cbs = [cb[i : i + cb_len] for i in range(0, len(cb), cb_len)]
    r1_seq = "".join([cbs[0], LINKER1, cbs[1], LINKER2, cbs[2], "C", umi, "T" * 18])
    r1_qual = "F" * len(r1_seq)
    r1 = fastq_str(query_name, r1_seq, r1_qual)

    r2_seq = segment.get_forward_sequence()
    r2_qual = "".join(map(lambda x: chr(x + 33), segment.get_forward_qualities()))
    r2 = fastq_str(query_name, r2_seq, r2_qual)
    return r1, r2


def get_cb_len(bam_file):
    with pysam.AlignmentFile(bam_file, "rb") as f:
        for segment in f:
            attr = segment.query_name.split("_")
            if len(attr) != 3:
                raise ValueError(
                    f"Currently only support celescope V1 BAM. {segment.query_name} is not a valid query name"
                )
            cb = attr[0]
            if len(cb) not in (24, 27):
                raise ValueError(f"{cb} not a valid cb")
            return len(cb) // 3


def main():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bam", required=True)
    parser.add_argument("-s", "--sample")
    parser.add_argument("-o", "--outdir", default="./")
    args = parser.parse_args()

    sample = args.sample if args.sample else os.path.basename(args.bam).split("_")[0]
    cb_len = get_cb_len(args.bam)
    f1_fn = f"{args.outdir}/{sample}_R1.fastq.gz"
    f2_fn = f"{args.outdir}/{sample}_R2.fastq.gz"
    n = 0
    mod = 1000000
    with pysam.AlignmentFile(args.bam, "rb") as f:
        with gzip.open(f1_fn, "wt") as f1, gzip.open(f2_fn, "wt") as f2:
            print("writing fastq...")
            for segment in f:
                r1, r2 = seg2fastq(segment, cb_len)
                f1.write(r1)
                f2.write(r2)
                n += 1
                if n % mod == 0:
                    print(f"{n//mod}M reads processed")


if __name__ == "__main__":
    main()
