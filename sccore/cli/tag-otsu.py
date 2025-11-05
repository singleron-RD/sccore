import glob
import os
from pathlib import Path
from celescope.tools.capture.threshold import Threshold
from celescope.tools.matrix import CountMatrix
import pandas as pd
import multiprocessing as mp


def run_single(celescope_dir, sample):
    """
    处理单个样本的函数
    """
    outdir = Path(sample)
    outdir.mkdir(exist_ok=True)
    mtx_path = f"{celescope_dir}/{sample}_c3/{sample}/03.count_cite/{sample}_citeseq.mtx.gz"
    df = pd.read_csv(mtx_path, sep="\t", index_col=0)
    array = df.iloc[0]
    obj = Threshold(array, threshold_method="otsu", otsu_plot_path=outdir / f"{sample}_tag_otsu.png", log_base=2)
    t = obj.run()
    pos_barcodes = array[array > t].index.tolist()
    neg_barcodes = array[array <= t].index.tolist()

    # write matrix
    mapfile_path = glob.glob(f"{celescope_dir}/{sample}_c3/*.mapfile")[0]
    match_dir = pd.read_csv(mapfile_path, sep="\t", index_col=False, header=None).iloc[0, 3]
    mtx_path = f"{match_dir}/outs/filtered"
    mtx = CountMatrix.from_matrix_dir(mtx_path)
    pos_mtx = mtx.slice_matrix_bc(pos_barcodes)
    neg_mtx = mtx.slice_matrix_bc(neg_barcodes)
    pos_mtx.to_matrix_dir(outdir / f"{sample}_tag_positive_matrix")
    neg_mtx.to_matrix_dir(outdir / f"{sample}_tag_negative_matrix")

    return t, len(array), len(pos_barcodes), len(neg_barcodes)


def run_single_wrapper(args):
    """
    包装函数，用于多进程映射
    """
    celescope_dir, sample = args
    return sample, run_single(celescope_dir, sample)


def find_samples_in_directory(celescope_dir):
    """
    在单个目录中查找所有样本
    """
    samples = []
    for d in os.listdir(celescope_dir):
        if d.endswith("_c3"):
            sample_name = d.replace("_c3", "")
            samples.append((celescope_dir, sample_name))
    return samples


def run_all_samples(celescope_dirs, num_processes=None):
    """
    使用多进程处理多个目录中的所有样本

    参数:
    celescope_dirs: 字符串或字符串列表，包含Celescope数据的目录路径
    num_processes: 进程数，默认为CPU核心数的一半
    """
    # 统一处理为列表
    if isinstance(celescope_dirs, str):
        celescope_dirs = [celescope_dirs]

    # 在所有目录中查找样本
    all_samples = []
    for celescope_dir in celescope_dirs:
        print(f"[INFO] Searching for samples in: {celescope_dir}")
        samples_in_dir = find_samples_in_directory(celescope_dir)
        all_samples.extend(samples_in_dir)
        print(f"[INFO] Found {len(samples_in_dir)} samples in {celescope_dir}")

    if not all_samples:
        print("[WARNING] No samples found in any directory!")
        return pd.DataFrame()

    print(f"[INFO] Total samples to process: {len(all_samples)}")

    # 设置进程数
    if num_processes is None:
        num_processes = max(1, mp.cpu_count() // 2)

    print(f"[INFO] Using {num_processes} processes")

    # 使用进程池
    with mp.Pool(processes=num_processes) as pool:
        results = {}
        # 使用imap_unordered获取完成的结果
        for sample, result in pool.imap_unordered(run_single_wrapper, all_samples):
            t, total_cells, pos_cells, neg_cells = result
            results[sample] = {
                "Sample": sample,
                "Total Cells": total_cells,
                "Positive Cells": pos_cells,
                "Negative Cells": neg_cells,
                "OTSU Threshold": t,
            }
            print(f"[PROGRESS] Completed {sample}: {pos_cells}/{total_cells} positive cells")

    # 转换为DataFrame
    records = list(results.values())
    df_summary = pd.DataFrame(records)
    return df_summary


if __name__ == "__main__":
    # 多个数据目录
    celescope_dirs = [
        "/SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/P25090402_cite/20251024",
        "/SGRNJ07/Standard_Analysis/PROJ03/PROJ_23.Other/P25090402_cite/20251027",
    ]

    # 使用多进程版本
    print("[INFO] Starting multiprocessing analysis...")
    df_summary = run_all_samples(celescope_dirs, num_processes=16)

    # 保存结果
    df_summary.to_csv("celescope_tag_summary.tsv", sep="\t", index=False)
    print("[INFO] Analysis completed!")
    print(df_summary)
