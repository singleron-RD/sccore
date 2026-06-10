#!/usr/bin/env python3
"""
create_mirna_reference.py

从 miRBase 的 mature.fa 文件中提取指定物种的 miRNA 序列，
生成 STARsolo 所需的 FASTA 和 GTF 参考文件。

用法示例:
    python create_mirna_reference.py -s hsa -i mature.fa -o mirna_ref

输入文件格式 (miRBase mature.fa):
    >hsa-let-7a-5p MIMAT0000062 Homo sapiens let-7a-5p
    UGAGGUAGUAGGUUGUAUAGUU
    >hsa-let-7a-3p MIMAT0004481 Homo sapiens let-7a-3p
    CUAUACAAUCUACUGUCUUUC
    ...
"""

import os
import argparse


def create_mirna_reference(species_code, input_fasta, output_dir="mirna_ref"):
    """
    从 miRBase mature.fa 创建适用于 STARsolo 的 fasta 和 gtf 文件。

    参数:
        species_code (str): 物种缩写 (如 'cel', 'hsa')
        input_fasta (str): miRBase 下载的 mature.fa 文件路径
        output_dir (str): 输出文件夹路径 (默认: mirna_ref)

    返回:
        tuple: (输出FASTA路径, 输出GTF路径)
    """
    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    out_fa_path = os.path.join(output_dir, f"{species_code}_mature.fasta")
    out_gtf_path = os.path.join(output_dir, f"{species_code}_mature.gtf")

    # 用于记录是否属于目标物种以及积累序列
    is_target = False
    current_mirna_id = ""  # 例如: hsa-let-7a-5p
    current_mimat_id = ""  # 例如: MIMAT0000062
    current_seq_lines = []  # 存储序列行

    with open(input_fasta, "r") as f_in, open(out_fa_path, "w") as f_out, open(out_gtf_path, "w") as gtf_out:
        for line in f_in:
            line = line.rstrip("\n")
            if line.startswith(">"):
                # 遇到新的 miRNA 条目：先处理前一条（如果有）
                if is_target and current_seq_lines:
                    process_mirna(current_mirna_id, current_mimat_id, current_seq_lines, f_out, gtf_out)

                # 解析头部
                header = line[1:].strip()  # 去掉 '>'
                parts = header.split()
                if len(parts) < 2:
                    # 格式异常，跳过
                    is_target = False
                    continue

                mirna_id = parts[0]  # 例如 hsa-let-7a-5p
                mimat_id = parts[1]  # 例如 MIMAT0000062

                # 检查是否为目标物种
                if mirna_id.startswith(f"{species_code}-"):
                    is_target = True
                    current_mirna_id = mirna_id
                    current_mimat_id = mimat_id
                    current_seq_lines = []
                else:
                    is_target = False
                    current_seq_lines = []
            else:
                # 非空行视为序列行（跳过空行）
                if line.strip() and is_target:
                    # 将 U 转换为 T（DNA 格式）
                    seq_line = line.strip().replace("U", "T").replace("u", "t")
                    current_seq_lines.append(seq_line)

        # 处理最后一条 miRNA
        if is_target and current_seq_lines:
            process_mirna(current_mirna_id, current_mimat_id, current_seq_lines, f_out, gtf_out)

    print(f"处理完成！\nFASTA: {out_fa_path}\nGTF:   {out_gtf_path}")
    return out_fa_path, out_gtf_path


def process_mirna(mirna_id, mimat_id, seq_lines, fa_handle, gtf_handle):
    """
    将单条 miRNA 写入 FASTA 和 GTF 文件。

    参数:
        mirna_id (str): miRNA 名称 (如 hsa-let-7a-5p)
        mimat_id (str): MIMAT 登录号
        seq_lines (list): 序列行列表
        fa_handle: FASTA 文件句柄
        gtf_handle: GTF 文件句柄
    """
    full_seq = "".join(seq_lines)
    seq_len = len(full_seq)
    if seq_len == 0:
        return

    # 1. 写入 FASTA
    fa_handle.write(f">{mirna_id}\n{full_seq}\n")

    # 2. 构建 GTF 属性字符串
    #    注意：属性值用双引号包裹，末尾有分号
    attr = (
        f'gene_id "{mimat_id}"; ' f'gene_name "{mirna_id}"; ' f'transcript_id "{mimat_id}.1"; ' f'gene_biotype "miRNA";'
    )

    # GTF 的起始、结束坐标均为 1 到序列长度（整条序列视作一个外显子）
    start = 1
    end = seq_len
    strand = "+"  # miRNA 通常不加链向信息，STARsolo 默认正链即可
    source = "miRBase"
    frame = "."

    # 写入 gene 行
    gtf_handle.write(f"{mirna_id}\t{source}\tgene\t{start}\t{end}\t{frame}\t{strand}\t{frame}\t{attr}\n")
    # 写入 transcript 行
    gtf_handle.write(f"{mirna_id}\t{source}\ttranscript\t{start}\t{end}\t{frame}\t{strand}\t{frame}\t{attr}\n")
    # 写入 exon 行
    gtf_handle.write(f"{mirna_id}\t{source}\texon\t{start}\t{end}\t{frame}\t{strand}\t{frame}\t{attr}\n")


def main():
    parser = argparse.ArgumentParser(description="从 miRBase mature.fa 生成 STARsolo 用的 miRNA 参考文件")
    parser.add_argument("-s", "--species", required=True, help="物种缩写 (如 'hsa', 'cel', 'mmu')")
    parser.add_argument("-i", "--input", required=True, help="miRBase 的 mature.fa 文件路径")
    parser.add_argument("-o", "--output-dir", default="mirna_ref", help="输出目录 (默认: mirna_ref)")
    args = parser.parse_args()

    # 检查输入文件是否存在
    if not os.path.isfile(args.input):
        print(f"错误：输入文件不存在 - {args.input}")
        return 1

    # 调用主函数
    create_mirna_reference(args.species, args.input, args.output_dir)
    return 0


if __name__ == "__main__":
    exit(main())
