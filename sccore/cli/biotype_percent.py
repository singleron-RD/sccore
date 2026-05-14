"""
输出两个tsv统计文件
1. 每个基因的gene_biotype, UMI读数和百分比，从高到低排序
2. 每种gene_biotype的总计UMI读数和百分比，从高到低排序
"""

import numpy as np

import scanpy as sc
import pandas as pd
import argparse
import os
from glob import glob


def parse_gtf(gtf_file):
    gene_info = {}
    with open(gtf_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "gene":
                continue
            info_field = fields[8]
            info_dict = {}
            for item in info_field.split(";"):
                item = item.strip()
                if item == "":
                    continue
                key, val = item.split(" ")
                info_dict[key] = val.strip('"')
            gene_name = info_dict.get("gene_name")
            if not gene_name:
                gene_name = info_dict.get("gene_id", "NA")
            gene_biotype = info_dict.get("gene_biotype", "NA")
            gene_info[gene_name] = gene_biotype
    return gene_info


def process_matrix(matrix_dir, gene_info, out_prefix):
    print(f"Processing sample: {out_prefix}")
    adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", make_unique=True)

    counts = np.array(adata.X.sum(axis=0)).ravel()  # 不会转稠密，只返回一个小的一维数组
    gene_counts = pd.DataFrame(
        {
            "gene_symbol": adata.var_names,
            "UMI_count": counts,
        }
    )

    gene_counts["gene_biotype"] = gene_counts["gene_symbol"].map(gene_info).fillna("NA")

    total_umi = gene_counts["UMI_count"].sum()
    gene_counts["UMI_percent"] = gene_counts["UMI_count"] / total_umi * 100
    gene_counts = gene_counts.sort_values("UMI_count", ascending=False)
    gene_counts.to_csv(f"{out_prefix}_gene_stats.tsv", sep="\t", index=False)

    biotype_counts = gene_counts.groupby("gene_biotype")["UMI_count"].sum().reset_index()
    biotype_counts["UMI_percent"] = biotype_counts["UMI_count"] / biotype_counts["UMI_count"].sum() * 100
    biotype_counts = biotype_counts.sort_values("UMI_count", ascending=False)
    biotype_counts.to_csv(f"{out_prefix}_biotype_stats.tsv", sep="\t", index=False)

    print(f"Done for sample: {out_prefix}\n")


def main():
    parser = argparse.ArgumentParser(description="统计10X matrix的UMI按gene和gene_biotype")
    parser.add_argument("--matrix", help="10X matrix目录路径 (matrix.mtx + barcodes.tsv + genes.tsv/ features.tsv)")
    parser.add_argument("--gtf", required=True, help="基因注释GTF文件")
    parser.add_argument("--out", help="输出文件前缀（如果使用--matrix）")
    parser.add_argument("-c", "--celescope_dir", help="Celescope样本目录，自动扫描每个sample的outs/filtered")
    args = parser.parse_args()

    print(f"解析GTF文件{args.gtf}...")
    gene_info = parse_gtf(args.gtf)

    if args.celescope_dir:
        samples_dirs = glob(os.path.join(args.celescope_dir, "*", "outs", "filtered"))
        if not samples_dirs:
            print(f"未找到任何 filtered 目录在 {args.celescope_dir}")
            return
        for matrix_dir in samples_dirs:
            sample_name = os.path.basename(os.path.dirname(os.path.dirname(matrix_dir)))  # 获取 sample 名
            process_matrix(matrix_dir, gene_info, sample_name)
    else:
        if not args.matrix or not args.out:
            print("如果不使用 -c，需要同时提供 --matrix 和 --out")
            return
        process_matrix(args.matrix, gene_info, args.out)


if __name__ == "__main__":
    main()
