#!/usr/bin/env python

import argparse
import gzip
import random
from collections import defaultdict

import pysam
import pandas as pd


def openfile(file_name, mode="rt", **kwargs):
    """open gzip or plain file"""
    if file_name.endswith(".gz"):
        file_obj = gzip.open(file_name, mode=mode, **kwargs)
    else:
        file_obj = open(file_name, mode=mode, **kwargs)
    return file_obj


def read_one_col(fn):
    """read one column file into list"""
    with openfile(fn) as f:
        return [x.strip() for x in f]


def get_records(bam_file, barcodes):
    a = []
    n_read = 0
    dup_align_read_names = set()
    with pysam.AlignmentFile(bam_file) as bam:
        for record in bam:
            n_read += 1
            if n_read % 1000000 == 0:
                print(f"processed {n_read} reads")
            cb = record.get_tag("CB")
            ub = record.get_tag("UB")
            try:
                gx = record.get_tag("GX")
            except Exception:
                continue
            if cb not in barcodes:
                continue
            if all(x != "-" for x in (cb, ub, gx)):
                if record.get_tag("NH") > 1:
                    if record.query_name in dup_align_read_names:
                        continue
                    else:
                        dup_align_read_names.add(record.query_name)
                # use int instead of str to avoid memory hog
                a.append((cb, ub, gx))
    return a


def sub_matrix(reads):
    subsamples = [1.0, 0.5, 0.1]
    n_reads = len(reads)
    expr_data = defaultdict(lambda: defaultdict(int))

    for frac in subsamples:
        # Step 1: 直接取前 frac 部分 reads
        cutoff = int(n_reads * frac)
        sampled_reads = reads[:cutoff]

        # Step 2: UMI 去重
        umi_sets = defaultdict(set)  # key: (barcode, gene), value: set of UMIs
        for barcode, umi, gene in sampled_reads:
            umi_sets[(barcode, gene)].add(umi)

        # Step 3: 统计表达量
        for (barcode, gene), umi_set in umi_sets.items():
            expr_data[gene][f"{barcode}_sub{frac}"] = len(umi_set)

    # 转换成 DataFrame
    expr_matrix = pd.DataFrame.from_dict(expr_data, orient="index").fillna(0).astype(int)
    expr_matrix = expr_matrix.reindex(sorted(expr_matrix.columns), axis=1)
    return expr_matrix


def main(args):
    barcodes = set(read_one_col(args.cell_barcode))
    a = get_records(args.bam, barcodes)
    random.seed(0)
    random.shuffle(a)

    expr_matrix = sub_matrix(a)
    expr_matrix.to_csv(f"{args.sample}_sub_matrix.tsv", sep="\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="saturation")
    parser.add_argument("-b", "--bam", help="bam file", required=True)
    parser.add_argument("-c", "--cell_barcode", help="barcode file", required=True)
    parser.add_argument("-s", "--sample", help="sample name", required=True)
    args = parser.parse_args()
    main(args)
