#!usr/bin/env python3
"""
https://gitee.com/singleron-rd/datasets/raw/main/10k_mouse_bone_marrow_V1/mouse_bone_marrow_filtered_feature_bc_matrix.tar.gz

https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/Sample_Y_S1_L001_R1_001.fastq.gz
"""

import argparse
import subprocess
import os
import csv


def run(fq_dict, n, sample, outf, prefix, repeat=1, skip=1e5, dry_run=False, has_fq3=False):
    """
    skip some reads, then get the following n reads from fq{pair}, write to {sample}_00{repeat}_R{pair}.fq.gz
    """
    n = int(n)
    skip = int(skip)
    nums = [1, 2]
    if has_fq3:
        nums.append(3)
    for r in range(1, repeat + 1):
        fqs = []
        for i in nums:
            total = (skip + n * r) * 4
            fq_name = f"{sample}_00{r}_R{i}.fq.gz"
            fqs.append(prefix + fq_name)
            cmd = f"zcat {fq_dict[i]} | head -n {total} | tail -n {n*4} | gzip > {fq_name}"
            print(cmd)
            if not dry_run:
                subprocess.check_call(cmd, shell=True)
        csv_line = ",".join([sample] + fqs) + "\n"
        print(csv_line)
        outf.write(csv_line)


def main():
    parser = argparse.ArgumentParser()
    # add: fq1, fq2, github or gitee, repo_name
    parser.add_argument("--fq1", help="fastq1", required=True)
    parser.add_argument("--fq2", help="fastq2", required=True)
    parser.add_argument("--fq3", help="fastq3")
    parser.add_argument("--website", "-w", help="github or gitee", required=True)
    parser.add_argument("--repo_name", "-r", help="repo name", required=True)
    parser.add_argument("--dir_name", "-d", help="dir name", required=True)
    parser.add_argument("--n_reads", "-n", help="number of reads in each fastq", default=5 * 1e4, type=int)
    # add dry_run
    parser.add_argument("--dry_run", "-dr", help="dry run", action="store_true")
    args = parser.parse_args()

    if args.website == "github":
        org = "singleron-RD"
    else:
        org = "singleron-rd"
    prefix = f"https://{args.website}.com/{org}/{args.repo_name}/raw/master/{args.dir_name}/"
    fq_dict = {1: args.fq1, 2: args.fq2}
    if not os.path.exists(args.dir_name):
        os.mkdir(args.dir_name)
    os.chdir(args.dir_name)
    outf = open("samplesheet.csv", "w")
    colnames = ["sample", "fastq_1", "fastq_2"]
    has_fq3 = False
    if args.fq3:
        fq_dict[3] = args.fq3
        colnames.append("fastq_3")
        has_fq3 = True
    outf.write(",".join(colnames) + "\n")
    run(fq_dict, args.n_reads, "sampleX", outf, prefix, repeat=2, dry_run=args.dry_run, has_fq3=has_fq3)
    run(fq_dict, args.n_reads, "sampleY", outf, prefix, repeat=1, dry_run=args.dry_run, has_fq3=has_fq3)
    outf.close()

    manifest_file = "manifest.csv"
    with open(manifest_file, "w") as f:
        fh = csv.DictWriter(f, fieldnames=["sample", "prefix"])
        fh.writeheader()
        fh.writerow({"sample": "sampleX", "prefix": "sampleX"})
        fh.writerow({"sample": "sampleY", "prefix": "sampleY"})


if __name__ == "__main__":
    main()
