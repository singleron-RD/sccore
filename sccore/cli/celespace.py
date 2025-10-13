import argparse
import glob
from pathlib import Path
import scanpy as sc
import pandas as pd
import numpy as np
import os
import matplotlib.image as mpimg
import json
import matplotlib.pyplot as plt
import matplotlib.colors as colors

from sccore import utils

SPATIAL_DIR = Path("/SGRNJ06/randd/USER/zhouyiqi/work/analysis/space/fake_spatial/")


@utils.add_log
def run_scanpy(sample, mtx_path, spatial_path):
    sample_path = Path(sample)
    if not sample_path.exists():
        os.makedirs(sample_path, exist_ok=True)
    print(f"Processing sample: {sample}\n")
    adata = sc.read_10x_mtx(mtx_path, var_names="gene_symbols", make_unique=True)

    pos_file = spatial_path / "tissue_positions_list.csv"
    positions = pd.read_csv(pos_file, header=None)
    positions.columns = ["barcode", "in_tissue", "array_row", "array_col", "pxl_col_in_fullres", "pxl_row_in_fullres"]
    # 保证 barcode 顺序和表达矩阵一致
    positions = positions.set_index("barcode").loc[adata.obs_names]
    adata.obs = pd.concat([adata.obs, positions], axis=1)
    adata.obsm["spatial"] = positions[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()

    sf_file = spatial_path / "scalefactors_json.json"
    with open(sf_file, "r") as f:
        scalefactors = json.load(f)
    hires = mpimg.imread(spatial_path / "tissue_hires_image.png")
    hires = np.stack([hires] * 3, axis=-1)
    lowres = mpimg.imread(spatial_path / "tissue_lowres_image.png")
    lowres = np.stack([lowres] * 3, axis=-1)
    adata.uns["spatial"] = {"sample1": {"images": {"hires": hires, "lowres": lowres}, "scalefactors": scalefactors}}

    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

    # 1. QC 图片
    plt.rcParams["figure.figsize"] = (8, 8)
    sc.pl.spatial(
        adata,
        img_key="hires",
        color=["total_counts", "n_genes_by_counts"],
        color_map="jet",
        size=1.5,
        norm=colors.LogNorm(vmin=1),
        show=False,
        save=None,
    )
    plt.savefig(sample_path / f"{sample}_qc_metrics.png", dpi=300, bbox_inches="tight")
    plt.close()

    # 2. cluster 图片
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.8)

    plt.rcParams["figure.figsize"] = (8, 8)
    sc.pl.spatial(
        adata,
        color=["leiden"],
        img_key="hires",
        size=1.5,
        show=False,
        cmap="tab20",
        save=None,
    )
    plt.savefig(sample_path / f"{sample}_cluster.png", dpi=300, bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--celescope_dir", type=str, required=True)
    parser.add_argument("-b", "--barcode_dim", type=int, required=True, choices=[24, 96])
    args = parser.parse_args()

    spatial_path = SPATIAL_DIR / f"barcode{args.barcode_dim}-spatial"
    for d in args.celescope_dir.split(","):
        print(f"checking path: {d}\n")
        mtx_path = glob.glob(f"{d}/*/*/raw")
        for path in mtx_path:
            sample = path.strip(d).split("/")[0]
            run_scanpy(sample, path, spatial_path)


if __name__ == "__main__":
    main()
