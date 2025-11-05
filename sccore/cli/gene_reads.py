#!/usr/bin/env python3
from pathlib import Path
import pysam
import pandas as pd
from collections import defaultdict
import sys


def count_gene_reads_and_umis(bam_path) -> pd.DataFrame:
    """
    统计 unique-mapped reads (NH:i:1) 的每个基因的 read 数量 和 unique UMI 数量（考虑 cell barcode）。
    返回 pandas DataFrame。
    """
    bam = pysam.AlignmentFile(bam_path, "rb")

    gene_reads = defaultdict(int)
    gene_umis = defaultdict(set)  # (cell, umi) 去重

    for read in bam:
        if read.is_unmapped:
            continue

        # unique mapping 判定
        try:
            if read.get_tag("NH") != 1:
                continue
        except KeyError:
            continue

        # 获取标签
        tags = dict(read.get_tags())
        gene = tags.get("GN")
        umi = tags.get("UB")
        cell = tags.get("CB")

        if not gene or not umi or not cell or gene == "-" or umi == "-" or cell == "-":
            continue

        gene_reads[gene] += 1
        gene_umis[gene].add((cell, umi))

    bam.close()

    # 生成 DataFrame
    df = pd.DataFrame(
        [(g, gene_reads[g], len(gene_umis[g])) for g in gene_reads], columns=["gene", "read_count", "umi_count"]
    )

    # 计算总数并添加百分比列（保留三位小数）
    total_reads = df["read_count"].sum()
    total_umis = df["umi_count"].sum()

    df["read_percent"] = (df["read_count"] / total_reads * 100).round(3)
    df["umi_percent"] = (df["umi_count"] / total_umis * 100).round(3)

    # 默认按 read_count 降序排列
    df = df.sort_values(by="read_count", ascending=False, ignore_index=True)
    return df


def main():
    if len(sys.argv) < 2:
        print("Usage: python count_gene_reads_umis_df.py <input.bam>")
        sys.exit(1)
    bam_path = Path(sys.argv[1])
    prefix = bam_path.stem
    out_path = Path(f"{prefix}.gene_read_umi.tsv")

    print("🔍 开始统计基因的 reads 和 UMIs...")
    df = count_gene_reads_and_umis(bam_path)
    df.to_csv(out_path, sep="\t", index=False)

    print(f"✅ 输出完成: {out_path} ({len(df)} genes)")
    print(df.head(10))  # 显示前10个结果


if __name__ == "__main__":
    main()
