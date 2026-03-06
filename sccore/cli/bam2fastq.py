from pathlib import Path
import sys
import pysam


def extract_unmapped_to_fastq(bam_path, r1_path, r2_path):
    # 打开 BAM 文件
    bam = pysam.AlignmentFile(bam_path, "rb")

    # 打开输出文件
    with open(r1_path, "w") as r1, open(r2_path, "w") as r2:
        for read in bam:
            # 仅处理未比对的 read
            if read.is_unmapped:
                # 获取标签：CB (Cell Barcode), UB (UMI)
                # 注意：有些流程中 UMI 标签可能是 'UR'
                raw_cb = read.get_tag("CB")
                if raw_cb == "-":
                    continue  # 跳过无效的 Barcode
                umi = read.get_tag("UB")

                # 处理 Barcode：去除下划线 (例如 "ATCG_1" -> "ATCG1")
                clean_cb = raw_cb.replace("_", "")

                # 构建 R1 的序列和质量值
                # R1 结构通常是: Barcode + UMI
                r1_seq = clean_cb + umi
                # 使用 'J' (Phred+33 质量值 41) 填充 R1 质量，或者根据需要截取
                r1_qual = "J" * len(r1_seq)

                # 获取 R2 的信息 (原始测序序列)
                r2_seq = read.query_sequence
                r2_qual = pysam.qualities_to_qualitystring(read.query_qualities)

                # 写入 R1 FASTQ
                r1.write(f"@{read.query_name}\n{r1_seq}\n+\n{r1_qual}\n")

                # 写入 R2 FASTQ
                r2.write(f"@{read.query_name}\n{r2_seq}\n+\n{r2_qual}\n")

    bam.close()
    print(f"提取完成！R1: {r1_path}, R2: {r2_path}")


if __name__ == "__main__":
    bam_file = sys.argv[1]
    sample = Path(bam_file).stem.split("_")[0]
    r1_file = f"{sample}_R1.fastq"
    r2_file = f"{sample}_R2.fastq"
    extract_unmapped_to_fastq(bam_file, r1_file, r2_file)
