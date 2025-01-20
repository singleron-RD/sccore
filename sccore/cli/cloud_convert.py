import os
import tarfile
import glob


def extract_tar_files(dir):
    # 在路径目录下查找符合 "_EmptyDrops_CR_matrix_10X.tar" 的文件
    tar_files = glob.glob(os.path.join(dir, "call-starsolo/execution/*_EmptyDrops_CR_matrix_10X.tar"))

    for tar_file in tar_files:
        # 从文件路径中提取文件名
        filename = os.path.basename(tar_file)
        # 提取 sample 名称
        sample = filename.split("_EmptyDrops_CR_matrix_10X.tar")[0]
        # 构造解压目录路径
        extract_path = os.path.join(sample, "outs", "filtered")
        # 确保解压目录存在
        os.makedirs(extract_path, exist_ok=True)

        # 解压 tar 文件
        with tarfile.open(tar_file, "r") as tar:
            tar.extractall(path=extract_path)
        print(f"Extracted {filename} to {extract_path}")


def copy_analysis_files(dir):
    tsne_file = glob.glob(os.path.join(dir, "call-analysis/execution/*_EmptyDrops_CR_tsne_coord.tsv"))[0]
    markers_file = glob.glob(os.path.join(dir, "call-analysis/execution/*_EmptyDrops_CR_markers.tsv"))[0]
    sample = (os.path.basename(tsne_file)).split("_EmptyDrops_CR_tsne_coord.tsv")[0]

    outdir = os.path.join(sample, "outs")
    os.makedirs(outdir, exist_ok=True)

    # copy to outdir
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
        if path:
            extract_tar_files(path)
            copy_analysis_files(path)
