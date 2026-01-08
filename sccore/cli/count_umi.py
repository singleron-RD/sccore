"""
输入一个fastq文件， fastq文件名称为{sample}_2.fq
read name格式为barcode:umi:id
统计每个barcode的read数目和umi数目，输出为{sample}_count.tsv文件，包含barcode, read_count, umi_count三列
"""

#!/usr/bin/env python3
import sys
import os

import pandas as pd
import pysam
from collections import defaultdict


def main():
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: python count_barcode_umi.py <input.fq>\n")
        sys.exit(1)

    fq = sys.argv[1]
    sample = os.path.basename(fq).split("_")[0]
    out_tsv = f"{sample}_count.tsv"

    read_count = defaultdict(int)
    umi_set = defaultdict(set)

    print(f"Processing file: {fq}")
    with pysam.FastxFile(fq) as f:
        for read in f:
            # read name format: barcode:umi:id
            try:
                barcode, umi, _ = read.name.split(":", 2)
            except ValueError:
                sys.stderr.write(f"Invalid read name format: {read.name}\n")
                continue

            read_count[barcode] += 1
            umi_set[barcode].add(umi)

    rows = []
    for barcode in read_count:
        rows.append({"barcode": barcode, "read_count": read_count[barcode], "umi_count": len(umi_set[barcode])})

    df = pd.DataFrame(rows)

    # 按 read_count 降序排序
    df = df.sort_values(by="read_count", ascending=False)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"Output written to: {out_tsv}")


if __name__ == "__main__":
    main()
