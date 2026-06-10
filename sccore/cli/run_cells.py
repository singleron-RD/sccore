#!/usr/bin/env python3
import argparse
import os
import subprocess
import concurrent.futures
import shutil


def find_samples(root_dir):
    """递归查找包含 outs 子目录的样本目录"""
    samples = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        if "outs" in dirnames:
            samples.append(dirpath)
    return samples


def run_cells(sample_dir, solo_cell_filter):
    """为单个样本运行 celescope rna cells"""

    sample_name = os.path.basename(sample_dir.rstrip("/"))
    parent_dir = os.path.dirname(sample_dir.rstrip("/"))
    outdir = os.path.join(sample_dir, "cells")

    if os.path.exists(outdir):
        print(f"[{sample_name}] cells directory already exists. removing: {outdir}", flush=True)
        shutil.rmtree(outdir)

    cmd = [
        "celescope",
        "rna",
        "cells",
        "--soloCellFilter",
        solo_cell_filter,
        "--outdir",
        f"{sample_name}/cells",
        "--sample",
        sample_name,
    ]

    print(f"[{sample_name}] Running celescope rna cells", flush=True)
    print(f"[{sample_name}] cd {parent_dir} && {' '.join(cmd)}", flush=True)

    _result = subprocess.run(
        cmd,
        cwd=parent_dir,
        check=True,
        text=True,
    )

    print(f"[{sample_name}] Completed", flush=True)
    return sample_name, "success"


def main():
    parser = argparse.ArgumentParser(description="Run celescope rna cells for multiple samples")
    parser.add_argument("root_dir", help="Root directory to search for samples")
    parser.add_argument(
        "--soloCellFilter",
        default="EmptyDrops_CR 3000 0.99 10 45000 90000 1000 0.01 20000 0.001 10000",
        help="soloCellFilter parameter for celescope rna cells (default: %(default)s)",
    )
    args = parser.parse_args()

    samples = find_samples(args.root_dir)

    if not samples:
        print("No samples found")
        return

    print(f"Found {len(samples)} samples:")
    for sample in samples:
        print(f"  - {os.path.basename(sample)}")

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [executor.submit(run_cells, sample, args.soloCellFilter) for sample in samples]

        for future in concurrent.futures.as_completed(futures):
            sample_name, status = future.result()
            print(f"[SUMMARY] {sample_name}: {status}", flush=True)


if __name__ == "__main__":
    main()
