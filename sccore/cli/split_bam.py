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


def split_bam_by_barcode(input_bam, barcode_file, output_dir, output_format):
    """
    按照 CB 标签拆分 BAM 文件，并支持 BAM 或 GZIP 压缩的 FASTQ 输出
    """
    # 加载 barcodes
    valid_barcodes = load_barcodes(barcode_file)
    print(f"Loaded {len(valid_barcodes)} barcodes from {barcode_file}")

    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 打开输入 BAM 文件
    bamfile = pysam.AlignmentFile(input_bam, "rb")

    # 字典用于存储每个 barcode 对应的 BAM/FASTQ 输出文件
    barcode_files = {}

    try:
        # 遍历每条 BAM 记录
        for read in bamfile:
            # 获取 CB 标签
            cb_tag = dict(read.get_tags()).get("CB")
            if cb_tag and cb_tag in valid_barcodes:
                if cb_tag not in barcode_files:
                    if output_format == "bam":
                        output_path = os.path.join(output_dir, f"{cb_tag}.bam")
                        barcode_files[cb_tag] = pysam.AlignmentFile(output_path, "wb", header=bamfile.header)
                    elif output_format == "fastq":
                        output_path = os.path.join(output_dir, f"{cb_tag}.fastq.gz")
                        barcode_files[cb_tag] = gzip.open(output_path, "wt")  # 以文本模式写入 gzip 压缩文件

                # 写入数据
                if output_format == "bam":
                    barcode_files[cb_tag].write(read)
                elif output_format == "fastq":
                    fastq_str = f"@{read.query_name}\n{read.query_sequence}\n+\n{''.join(chr(q + 33) for q in read.query_qualities)}\n"
                    barcode_files[cb_tag].write(fastq_str)

    finally:
        for cb_file in barcode_files.values():
            cb_file.close()
        bamfile.close()


def main():
    parser = argparse.ArgumentParser(description="Split BAM file by CB tag using a list of barcodes")
    parser.add_argument("-i", "--input", required=True, help="Input BAM file")
    parser.add_argument("-b", "--barcode", required=True, help="Barcode file (gzip compressed, one barcode per line)")
    parser.add_argument("-o", "--output", required=True, help="Output directory for split files")
    parser.add_argument(
        "-f", "--format", choices=["bam", "fastq"], default="bam", help="Output format: bam or fastq.gz"
    )

    args = parser.parse_args()

    split_bam_by_barcode(args.input, args.barcode, args.output, args.format)


if __name__ == "__main__":
    main()
