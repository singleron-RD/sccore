#!/usr/bin/env python3
import argparse
import os
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np


def get_top10_mirnas(adata):
    """从 AnnData 中提取 top10 microRNA 基因"""
    mirna_genes = [gene for gene in adata.var_names if gene.startswith("hsa-") or gene.startswith("mmu-")]
    if not mirna_genes:
        return None
    mean_expression = adata[:, mirna_genes].layers["counts"].mean(axis=0)
    mean_expression = np.asarray(mean_expression).ravel()
    gene_means = list(zip(mirna_genes, mean_expression))
    gene_means.sort(key=lambda x: x[1], reverse=True)
    top10_mirnas = [gene for gene, _ in gene_means[:10]]
    return top10_mirnas


def get_cluster_col(adata):
    """尝试从 adata.obs 中获取常用的 cluster 列名"""
    for col in ["leiden", "louvain", "seurat_clusters", "cluster"]:
        if col in adata.obs.columns:
            return col
    return None


def find_samples(root_dir):
    """递归查找包含 outs/spatial 子目录的样本目录"""
    samples = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        spatial_path = os.path.join(dirpath, "outs", "spatial")
        if os.path.isdir(spatial_path):
            samples.append(dirpath)
    return samples


def plot_mirna_spatial(adata, top10_mirnas, sample_name, output_dir):
    """画出top10 microRNA的空间分布"""
    try:
        output_path = os.path.join(output_dir, f"{sample_name}_top10_mirna_spatial.png")

        sc.pl.spatial(adata, color=top10_mirnas, layer="normalized", cmap="Reds", size=1.5)
        # 保存，可调整 dpi 和 bbox_inches
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()  # 关闭当前图像，避免在交互环境中重复显示

        print(f"Saved spatial plot: {output_path}")

    except Exception as e:
        print(f"Error processing spatial plot for {sample_name}: {e}")


def plot_mirna_violin(adata, top10_mirnas, sample_name, output_dir):
    """画出top10 microRNA在每个cluster的vlnplot"""
    try:
        cluster_col = get_cluster_col(adata)
        if cluster_col is None:
            print(f"No cluster column found in {sample_name}, skipping violin plot")
            return

        output_path = os.path.join(output_dir, f"{sample_name}_top10_mirna_violin.png")

        n_genes = len(top10_mirnas)
        ncols = 4
        nrows = (n_genes + ncols - 1) // ncols  # 计算所需行数

        # 创建子图网格
        fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5 * ncols, 4 * nrows))
        axes = axes.flatten()  # 将多维轴对象展平为一维，方便迭代

        # 循环绘制每个基因的小提琴图
        for i, gene in enumerate(top10_mirnas):
            sc.pl.violin(adata, keys=gene, groupby="cluster", layer="normalized", rotation=90, ax=axes[i], show=False)

        # 隐藏多余的子图（如果基因数量不能填满网格）
        for j in range(i + 1, len(axes)):
            axes[j].axis("off")

        plt.tight_layout()
        plt.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close()

        print(f"Saved violin plot: {output_path}")

    except Exception as e:
        print(f"Error processing violin plot for {sample_name}: {e}")


def main():
    parser = argparse.ArgumentParser(description="Plot top10 microRNA spatial distribution from Visium data")
    parser.add_argument("root_dir", help="Root directory to search for sample directories")
    parser.add_argument("-o", "--output", default="mirna_spatial_plots", help="Output directory for plots")
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    samples = find_samples(args.root_dir)
    print(f"Found {len(samples)} sample directories")

    for sample in samples:
        print(f"Processing {sample}")
        try:
            adata = sc.read_h5ad(os.path.join(sample, "outs/rna.h5ad"))
        except Exception as e:
            print(f"Error reading h5ad for {sample}: {e}")
            continue

        top10_mirnas = get_top10_mirnas(adata)
        if top10_mirnas is None:
            print(f"No microRNA genes found in {sample}")
            continue

        sample_name = os.path.basename(sample)
        print(f"{sample_name} Top 10 microRNA genes: {top10_mirnas}")

        plot_mirna_spatial(adata, top10_mirnas, sample_name, args.output)
        plot_mirna_violin(adata, top10_mirnas, sample_name, args.output)


if __name__ == "__main__":
    main()
