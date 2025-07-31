import argparse
import pysam
import pandas as pd
from celescope.tools import utils


def read_fastq_names_pysam(fq_path):
    read_names = set()

    with pysam.FastxFile(fq_path) as fastq_file:
        for record in fastq_file:
            name = record.name[:-5]
            read_names.add(name)

    return read_names


class BamUMIStats:
    def __init__(self, args):
        self.bam = args.bam
        self.bclist = args.bclist
        self.args = args
        self.set3 = set()
        self.set5 = set()
        self.total_count = 0
        self.count3 = 0
        self.count5 = 0
        self.intersec_count = 0

    @utils.add_log
    def __call__(self):
        p5_names = read_fastq_names_pysam(self.args.p5_fq)

        df = pd.read_csv(self.bclist, names=["bc"])
        cell_bc = set(df.bc)
        inbam = pysam.AlignmentFile(self.bam, "rb")

        for read in inbam:
            if read.is_secondary:
                continue
            cb = read.get_tag("CB")
            umi = read.get_tag("UB")
            gx = read.get_tag("GX")
            if cb == "-" or gx == "-":
                continue
            if cb not in cell_bc:
                continue
            self.total_count += 1
            if str(read.query_name) in p5_names:
                self.count5 += 1
                self.set5.add((cb, umi))
            else:
                self.count3 += 1
                self.set3.add((cb, umi))
            if self.total_count % 1000000 == 0:
                print(f"processed {self.total_count} reads")
                print(f"p3 reads: {len(self.set3)}")
                print(f"p5 reads: {len(self.set5)}")

        inbam.close()
        self.intersec_count = len(self.set3 & self.set5)

        with open("bc_umi_count.txt", "w") as fp:
            fp.write(f"total (cb,umi,gene) read count : {self.total_count}\n")
            fp.write(f"3p read count : {self.count3}\n")
            fp.write(f"5p read count : {self.count5}\n")
            fp.write(f"3p umi count : {len(self.set3)}\n")
            fp.write(f"5p umi count : {len(self.set5)}\n")
            fp.write(f"3p5p umi intersection : {self.intersec_count}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count UMI per CB by tag (3p/5p)")
    parser.add_argument("--bam", help="bam file", required=True)
    parser.add_argument("--bclist", help="cell barcode file", required=True)
    parser.add_argument("--p5_fq", required=True)
    args = parser.parse_args()
    BamUMIStats(args)()
