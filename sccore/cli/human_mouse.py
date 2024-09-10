"""
estimate ambient from human_mouse sample
"""

from sccore.matrix import CountMatrix
import pandas as pd
import argparse
import glob
import os
import functools


def get_metrics_dict(mtx_path, doublet_threshold):
    mtx = CountMatrix.from_matrix_dir(mtx_path)
    m = mtx.get_matrix()
    f = mtx.get_features()

    human = []
    mouse = []
    for i, gn in enumerate(f.gene_name):
        if gn == gn.upper():
            human.append(i)
        else:
            mouse.append(i)
    umi_h = m.tocsr()[human, :].sum(axis=0)
    umi_m = m.tocsr()[mouse, :].sum(axis=0)
    df_m = pd.DataFrame(umi_m.T, columns=["mouse"])
    df_h = pd.DataFrame(umi_h.T, columns=["human"])
    df = pd.merge(df_m, df_h, left_index=True, right_index=True)
    df["umi_sum"] = df.sum(axis=1)
    df["human_percent"] = df["human"] / df["umi_sum"] * 100
    df["mouse_percent"] = df["mouse"] / df["umi_sum"] * 100
    df["identity"] = "doublet"
    df.loc[df["human_percent"] < doublet_threshold, "identity"] = "mouse"
    df.loc[df["mouse_percent"] < doublet_threshold, "identity"] = "human"
    df.loc[df["identity"] == "human", "ambient_percent"] = df[df["identity"] == "human"].mouse_percent
    df.loc[df["identity"] == "mouse", "ambient_percent"] = df[df["identity"] == "mouse"].human_percent

    n_cell = df.shape[0]
    n_doublet = df[df["identity"] == "doublet"].shape[0]
    doublet_percent = round(n_doublet / n_cell * 100, 2)
    dict = {
        "n_cell": n_cell,
        "n_human_cell": df[df["identity"] == "human"].shape[0],
        "n_mouse_cell": df[df["identity"] == "mouse"].shape[0],
        "doublet_umi_percent_threshold": doublet_threshold,
        "n_doublet": n_doublet,
        "doublet_cell_percent": doublet_percent,
        "median_ambient_umi_percent": round(df["ambient_percent"].median(), 2),
        "mean_ambient_umi_percent": round(df["ambient_percent"].mean(), 2),
    }

    return dict


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--celescope_dir", type=str, required=True)
    parser.add_argument("--doublet_threshold", type=float, default=25.0)
    parser.add_argument("--out_prefix", type=str, default="human_mouse")
    args = parser.parse_args()

    cur_dir = os.getcwd()
    dfs = []
    for d in args.celescope_dir.split(","):
        os.chdir(d)
        paths = glob.glob("*/filtered")
        for path in paths:
            sample = path.split("/")[0]
            dict = get_metrics_dict(path, args.doublet_threshold)
            dfs.append(pd.DataFrame.from_dict(dict, orient="index", columns=[sample]))

    df = functools.reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), dfs)
    df = df.T
    df = df.astype(
        {
            "n_cell": int,
            "n_human_cell": int,
            "n_mouse_cell": int,
            "n_doublet": int,
        }
    )
    os.chdir(cur_dir)
    df.to_csv(f"{args.out_prefix}.tsv", sep="\t")


if __name__ == "__main__":
    main()
