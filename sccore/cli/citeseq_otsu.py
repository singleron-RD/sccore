#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd
from celescope.tools.capture.threshold import Otsu  # 假设你的 Otsu 类保存在 Otsu.py 里


def main():
    # 输入文件
    matrix_file = sys.argv[1]

    # 输出文件
    output_matrix_file = "filtered_test1_citeseq.tsv.gz"
    log_file = "otsu_filter_log.txt"

    # 读取 tsv.gz
    print("Reading matrix...")
    df = pd.read_csv(matrix_file, sep="\t", index_col=0, compression="gzip")

    print(f"Matrix shape: {df.shape} (features x cells)")

    # 打开 log 文件
    log_f = open(log_file, "w")
    log_f.write("Feature\tOTSU_threshold\tNum_cells_after_filtering\n")

    # 对每一行应用 OTSU
    print("Applying OTSU filtering...")
    for i, feature_name in enumerate(df.index):
        row_data = df.loc[feature_name].values.astype(float)
        otsu = Otsu(row_data, otsu_plot_path=f"{feature_name}_otsu.png")
        threshold = otsu.run()

        # 低于阈值的设为 0
        row_data[row_data < threshold] = 0
        num_cells_after = np.sum(row_data > 0)

        # 更新 DataFrame
        df.loc[feature_name] = row_data

        # 写 log
        log_f.write(f"{feature_name}\t{threshold}\t{num_cells_after}\n")

    log_f.close()

    # 保存过滤后的矩阵
    print("Saving filtered matrix...")
    df.to_csv(output_matrix_file, sep="\t", compression="gzip")

    print("Done.")


if __name__ == "__main__":
    main()
