#!/usr/bin/env python3


def revcomp(seq: str) -> str:
    """返回 DNA 序列的反向互补序列"""
    complement = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(complement)[::-1]


def process_file(input_file: str, output_file: str):
    with open(input_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            barcode = line.strip()
            if barcode:  # 跳过空行
                fout.write(revcomp(barcode) + "\n")


if __name__ == "__main__":
    import sys

    if len(sys.argv) != 3:
        print(f"用法: {sys.argv[0]} input.txt output.txt")
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]
    process_file(input_file, output_file)
