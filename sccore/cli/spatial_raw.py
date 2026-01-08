import glob
import sys

import shutil
import h5py
import numpy as np
from pathlib import Path
import scanpy as sc


def convert_10x_h5(mtx_dir, outfile, library_id="library0"):
    """
    Convert a 10X MTX folder to a Visium-compatible HDF5 file.
    Existing file will be overwritten.

    Args:
        mtx_dir (str or Path): Path to 10X-formatted mtx directory.
        outfile (str or Path): Output HDF5 file path.
        library_id (str): Library ID for Visium metadata.
    """
    mtx_dir = Path(mtx_dir)
    outfile = Path(outfile)
    if not outfile.suffix == ".h5":
        outfile = outfile.with_suffix(".h5")

    # 如果文件已存在，直接覆盖
    if outfile.exists():
        outfile.unlink()

    # 读取 MTX 文件
    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols", make_unique=True)

    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # 默认 feature_type 和 genome
    feature_types = (
        adata.var["feature_types"] if "feature_types" in adata.var.columns else ["Gene Expression"] * n_genes
    )
    genomes = adata.var["genome"] if "genome" in adata.var.columns else ["unknown"] * n_genes
    gene_ids = adata.var["gene_ids"] if "gene_ids" in adata.var.columns else adata.var_names

    # 创建 HDF5 文件
    with h5py.File(outfile, "w") as f:
        grp = f.create_group("matrix")

        # barcodes
        grp.create_dataset(
            "barcodes",
            data=np.array(adata.obs_names, dtype=f"|S{max(map(len, adata.obs_names))}"),
        )

        # counts
        grp.create_dataset("data", data=adata.X.data.astype(np.int32))
        grp.create_dataset("indices", data=adata.X.indices.astype(np.int32))
        grp.create_dataset("indptr", data=adata.X.indptr.astype(np.int32))
        grp.create_dataset("shape", data=np.array([n_genes, n_cells], dtype=np.int64))

        # features
        ftrs = grp.create_group("features")
        ftrs.create_dataset("id", data=np.array(gene_ids, dtype=f"|S{max(map(len, gene_ids))}"))
        ftrs.create_dataset(
            "name",
            data=np.array(adata.var_names, dtype=f"|S{max(map(len, adata.var_names))}"),
        )
        ftrs.create_dataset(
            "feature_type",
            data=np.array(feature_types, dtype=f"|S{max(map(len, feature_types))}"),
        )
        ftrs.create_dataset("genome", data=np.array(genomes, dtype=f"|S{max(map(len, genomes))}"))

        # Visium-specific metadata
        all_tag_keys = np.array([b"feature_type", b"genome", b"id", b"name", b"library_ids"])
        ftrs.create_dataset("_all_tag_keys", data=all_tag_keys)

        library_ids_dataset = np.array([library_id.encode()] * n_genes)
        ftrs.create_dataset("library_ids", data=library_ids_dataset)

        # HDF5 attributes required by scanpy.read_visium
        f.attrs["library_ids"] = np.array([library_id.encode()])
        f.attrs["chemistry_description"] = "Spatial3"
        f.attrs["software_version"] = "CustomPythonScript"


def main():
    spatial_dir = Path(sys.argv[1]) / "call-sample_starsolo_analysis/execution/"
    all_matrix = Path(glob.glob(str(spatial_dir / "*_all_matrix.tar"))[0])
    sample = all_matrix.stem.replace("_all_matrix", "")
    spatial_tar = spatial_dir / f"{sample}_spatial.tar"
    html = spatial_dir / f"{sample}_report.html"

    # untar files
    shutil.unpack_archive(all_matrix, "./raw")
    shutil.unpack_archive(spatial_tar, "./spatial")
    shutil.copy(html, f"./{sample}_report.html")


if __name__ == "__main__":
    main()
