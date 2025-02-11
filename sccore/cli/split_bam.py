import subprocess
import pysam
import gzip
import argparse
import os


def load_barcodes(barcode_file):
    """
    读取 barcode.tsv.gz 文件，返回一个包含所有 barcodes 的集合
    """
    barcodes = set()
    with gzip.open(barcode_file, "rt") as f:
        for line in f:
            barcodes.add(line.strip())
    return barcodes


def split_bam_by_barcode(input_bam, barcode_file, output_dir):
    """
    按照 CB 标签拆分 BAM 文件，只输出 barcode 文件中包含的条形码
    """
    # 加载 barcodes
    valid_barcodes = load_barcodes(barcode_file)
    print(f"Loaded {len(valid_barcodes)} barcodes from {barcode_file}")

    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 打开输入 BAM 文件
    bamfile = pysam.AlignmentFile(input_bam, "rb")

    # 字典用于存储每个 barcode 对应的 BAM 输出文件
    barcode_files = {}

    try:
        # 遍历每条 BAM 记录
        for read in bamfile:
            # 获取 CB 标签
            cb_tag = dict(read.get_tags()).get("CB")
            if cb_tag and cb_tag in valid_barcodes:
                # 如果该 barcode 尚未创建 BAM 文件，则创建
                if cb_tag not in barcode_files:
                    output_bam_path = os.path.join(output_dir, f"{cb_tag}.bam")
                    barcode_files[cb_tag] = pysam.AlignmentFile(output_bam_path, "wb", header=bamfile.header)

                # 将该记录写入对应的 BAM 文件
                barcode_files[cb_tag].write(read)

    finally:
        for cb_file in barcode_files.values():
            cb_file.close()
        bamfile.close()


MAX_OPEN_FILE = 50000


def set_ulimit():
    """
    尝试在脚本中设置 ulimit -n 的文件描述符限制
    """
    try:
        subprocess.run(f"ulimit -n {MAX_OPEN_FILE}", shell=True, check=True)
        print(f"File descriptor limit set to {MAX_OPEN_FILE}.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to set ulimit: {e}")


def main():
    parser = argparse.ArgumentParser(description="Split BAM file by CB tag using a list of barcodes")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-b", "--barcode", required=True, help="Barcode file (gzip compressed, one barcode per line)")
    parser.add_argument("-o", "--output", required=True, help="Output directory for split BAM files")

    args = parser.parse_args()

    set_ulimit()

    split_bam_by_barcode(args.input, args.barcode, args.output)


if __name__ == "__main__":
    main()
