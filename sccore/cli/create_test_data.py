#!usr/bin/env python3
"""
https://gitee.com/singleron-rd/datasets/raw/main/10k_mouse_bone_marrow_V1/mouse_bone_marrow_filtered_feature_bc_matrix.tar.gz

https://github.com/singleron-RD/scrna_test_data/raw/master/GEXSCOPE-V2/Sample_Y_S1_L001_R1_001.fastq.gz
"""
import argparse
import subprocess
import os

def run(fq_dict, n, sample, outf, prefix, repeat=1, skip = 1e5):
    """
    skip some reads, then get the following n reads from fq{pair}, write to {sample}_00{repeat}_R{pair}.fq.gz
    """
    n = int(n)
    skip = int(skip)
    for r in range(1, repeat+1):
        fqs = []
        for i in [1,2]:
            total = (skip + n * r) * 4
            fq_name = f"{sample}_00{r}_R{i}.fq.gz"
            fqs.append(prefix + fq_name)
            cmd = f"zcat {fq_dict[i]} | head -n {total} | tail -n {n*4} | gzip > {fq_name}"
            print(cmd)
            subprocess.check_call(cmd, shell=True)
        csv_line = ",".join([sample] + fqs) + '\n'
        print(csv_line)
        outf.write(csv_line)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # add: fq1, fq2, github or gitee, repo_name
    parser.add_argument("--fq1", help="fastq1", required=True)
    parser.add_argument("--fq2", help="fastq2", required=True)
    parser.add_argument("--website",'-w', help="github or gitee", required=True)
    parser.add_argument("--repo_name", '-r', help="repo name", required=True)
    parser.add_argument("--dir_name", '-d', help="dir name", required=True)
    parser.add_argument("--n_reads", '-n', help="number of reads in each fastq", default=5*1e4, type=int)
    args = parser.parse_args()

    if args.website == "github":
        org = 'singleron-RD'
    else:
        org = 'singleron-rd'
    prefix = f"https://{args.website}.com/{org}/{args.repo_name}/raw/master/{args.dir_name}/"
    fq_dict = {1:args.fq1, 2:args.fq2}
    if not os.path.exists(args.dir_name):
        os.mkdir(args.dir_name)
    os.chdir(args.dir_name)
    outf = open("samplesheet.csv", 'w')
    run(fq_dict, args.n_reads, 'sampleX', outf, prefix, repeat=2)
    run(fq_dict, args.n_reads, 'sampleY', outf, prefix, repeat=1)
    outf.close()