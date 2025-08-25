import pandas as pd
import subprocess
import argparse
import os


def write_mapfile():
    data = [
        ["X", "./fastqs", "sampleX"],
        ["Y", "./fastqs", "sampleY"],
    ]
    df = pd.DataFrame(data, columns=["prefix", "folder", "sample"])
    df.to_csv("mapfile.tsv", sep="\t", index=False, header=False)


def subset_fastq(in_file, start: int, n_reads: int, out_file):
    outdir = "fastqs"
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    out_file = os.path.join(outdir, out_file)
    start_line = start * 4
    read_line = n_reads * 4
    total_line = start_line + read_line
    cmd = f"zcat {in_file} | head -n {total_line} | tail -n {read_line} | gzip > {out_file}"
    print(cmd)
    subprocess.check_call(cmd, shell=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq1", required=True)
    parser.add_argument("--fq2", required=True)
    args = parser.parse_args()

    n_reads = 50000
    fq_dict = {
        "R1": args.fq1,
        "R2": args.fq2,
    }
    for pair, fq_file in fq_dict.items():
        subset_fastq(fq_file, 100000, n_reads, f"X_001_{pair}.fastq.gz")
        subset_fastq(fq_file, 100000 + n_reads, n_reads, f"X_002_{pair}.fastq.gz")
        subset_fastq(fq_file, 100000 + n_reads * 2, n_reads, f"Y_001_{pair}.fastq.gz")

    write_mapfile()


if __name__ == "__main__":
    main()
