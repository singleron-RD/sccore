import os
import tarfile
import glob


def extract_tar_files(dir):
    # 递归查找所有 *_EmptyDrops_CR_matrix_10X.tar 文件
    tar_files = glob.glob(os.path.join(dir, "**", "*EmptyDrops_CR_matrix_10X.tar"), recursive=True)
    for tar_file in tar_files:
        filename = os.path.basename(tar_file)
        sample = filename.replace("_EmptyDrops_CR_matrix_10X.tar", "")
        extract_path = os.path.join(sample, "outs", "filtered")
        os.makedirs(extract_path, exist_ok=True)

        with tarfile.open(tar_file, "r") as tar:
            tar.extractall(path=extract_path)
        print(f"Extracted {filename} to {extract_path}")


def copy_analysis_files(dir):
    # 递归查找 tsne 和 markers 文件
    tsne_files = glob.glob(os.path.join(dir, "**", "*_EmptyDrops_CR_tsne_coord.tsv"), recursive=True)
    markers_files = glob.glob(os.path.join(dir, "**", "*_EmptyDrops_CR_markers.tsv"), recursive=True)

    # 构建字典：sample -> {tsne, markers}
    sample_files = {}
    for tsne in tsne_files:
        sample = os.path.basename(tsne).split("_EmptyDrops_CR_tsne_coord.tsv")[0]
        sample_files.setdefault(sample, {})["tsne"] = tsne

    for markers in markers_files:
        sample = os.path.basename(markers).split("_EmptyDrops_CR_markers.tsv")[0]
        sample_files.setdefault(sample, {})["markers"] = markers

    for sample, files in sample_files.items():
        tsne_file = files.get("tsne")
        markers_file = files.get("markers")
        if tsne_file and markers_file:
            outdir = os.path.join(sample, "outs")
            os.makedirs(outdir, exist_ok=True)
            os.system(f"cp {tsne_file} {outdir}/tsne_coord.tsv")
            os.system(f"cp {markers_file} {outdir}/markers.tsv")
            print(f"Copied analysis files to {outdir}")


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Convert cloud results to local match_dir")
    parser.add_argument("dir_paths_txt", type=str, help="txt file contains directory to convert, one per line")
    args = parser.parse_args()

    with open(args.dir_paths_txt, "r") as f:
        paths = f.read().splitlines()

    for path in paths:
        path = path.strip()
        if os.path.exists(path):
            print(f"Converting {path}")
        else:
            print(f"{path} does not exist, skip")
        if path:
            extract_tar_files(path)
            copy_analysis_files(path)


if __name__ == "__main__":
    main()
