"""
将AccuraSCOPE RNA或者DNA的fastq拆分到孔的fastq
- 根据well name的mapping关系命名对应的fastq文件
- 输出指标：valid reads百分比，在孔中的reads百分比
"""

from collections import defaultdict
import pandas as pd
import pysam
import argparse
from sccore import utils, parse_chemistry
from pathlib import Path
import sys


def get_well_barcode(bc_file: str) -> dict[int, str]:
    barcodes = utils.one_col_to_list(bc_file)
    return {i: x for i, x in enumerate(barcodes, start=1)}


def get_barcode_name(
    bc_file: str,
    well_name_file: str,
) -> dict[str, str]:
    well_barcode = get_well_barcode(bc_file)
    well_name = utils.two_col_to_dict(well_name_file)
    barcode_name = {}
    for well, barcode in well_barcode.items():
        if well in well_name:
            barcode_name[barcode] = well_name[well]
        else:
            barcode_name[barcode] = f"noise_{well}"
    return barcode_name


class SplitRNA:
    def __init__(self, args):
        self.fq1 = args.fq1.split(",")
        self.fq2 = args.fq2.split(",")
        self.p3_bcumi_runner = parse_chemistry.BcUmi("AccuraSCOPE_RNA_3p")
        self.p5_bcumi_runner = parse_chemistry.BcUmi("AccuraSCOPE_RNA_5p", strict=True)
        _pattern_dict, p3_bc_file = parse_chemistry.get_pattern_dict_and_bc("AccuraSCOPE_RNA_3p")
        _pattern_dict, p5_bc_file = parse_chemistry.get_pattern_dict_and_bc("AccuraSCOPE_RNA_5p")
        self.p3_barcode_name = get_barcode_name(p3_bc_file[0], args.well_name)
        self.p5_barcode_name = get_barcode_name(p5_bc_file[0], args.well_name)
        self.auto_runner = parse_chemistry.AutoAccuraRNA(self.fq1)
        self.outdir = Path("./rna")

    @utils.add_log
    def run(self):
        out_handles = {}
        total_reads = 0
        p3_reads = 0
        p5_reads = 0
        signal_reads = 0
        name_read_count = defaultdict(lambda: defaultdict(int))
        namesheet_df_lines = []

        # mkdir
        fastq_dir = self.outdir / "fastqs"
        if not fastq_dir.exists():
            fastq_dir.mkdir(parents=True, exist_ok=True)

        for fq1, fq2 in zip(self.fq1, self.fq2):
            with pysam.FastxFile(fq1) as handle1, pysam.FastxFile(fq2) as handle2:
                for read1, read2 in zip(handle1, handle2):
                    total_reads += 1
                    if total_reads % 1_000_000 == 0:
                        sys.stderr.write(f"Processed {total_reads} reads...\n")
                    chemistry = self.auto_runner.seq_chemistry(read1.sequence)
                    if not chemistry:
                        continue
                    if chemistry == "AccuraSCOPE_RNA_3p":
                        valid, _corrected, corrected_seq, _umi = self.p3_bcumi_runner.get_bc_umi(read1.sequence)
                    else:
                        valid, _corrected, corrected_seq, _umi = self.p5_bcumi_runner.get_bc_umi(read1.sequence)
                    if not valid:
                        continue
                    barcode = corrected_seq
                    if chemistry == "AccuraSCOPE_RNA_3p":
                        p3_reads += 1
                        name = self.p3_barcode_name[barcode]
                        name_read_count[name]["p3"] += 1
                    else:
                        p5_reads += 1
                        name = self.p5_barcode_name[barcode]
                        name_read_count[name]["p5"] += 1

                    if not name.startswith("noise_"):
                        signal_reads += 1
                        if name not in out_handles:
                            out_file = fastq_dir / f"{name}.fq.gz"
                            namesheet_df_lines.append(
                                {
                                    "sample": name,
                                    "fastq_1": str(out_file.absolute()),
                                    "fastq_2": "",
                                    "strandedness": "auto",
                                }
                            )
                            out_handles[name] = utils.generic_open(out_file, "wt", compresslevel=1)

                        out_handles[name].write(utils.fastq_str(read2.name, read2.sequence, read2.quality))

        # output metrics
        metrics_file = self.outdir / "rna_metrics.txt"
        p3_reads_percent = round(p3_reads / total_reads * 100, 2) if total_reads > 0 else 0
        p5_reads_percent = round(p5_reads / total_reads * 100, 2) if total_reads > 0 else 0
        signal_reads_percent = round(signal_reads / total_reads * 100, 2) if total_reads > 0 else 0
        metrics_dict = {
            "total_reads": total_reads,
            "p3_reads": f"{p3_reads} ({p3_reads_percent}%)",
            "p5_reads": f"{p5_reads} ({p5_reads_percent}%)",
            "signal_reads": f"{signal_reads} ({signal_reads_percent}%)",
        }
        with open(metrics_file, "w") as mf:
            for k, v in metrics_dict.items():
                mf.write(f"{k}: {v}\n")

        # read count per name
        df = pd.DataFrame.from_dict(name_read_count, orient="index")
        df = df.reindex(columns=["p3", "p5"]).fillna(0).astype(int)
        df.sort_values(by=["p3", "p5"], ascending=False, inplace=True)
        name_count_file = self.outdir / "rna_read_count.txt"
        df.to_csv(name_count_file, index_label="name", sep="\t")

        # namesheet
        namesheet_file = self.outdir / "rna_samplesheet.csv"
        namesheet_df = pd.DataFrame(namesheet_df_lines)
        namesheet_df.to_csv(namesheet_file, index=False)

        # close file handles
        for handle in out_handles.values():
            handle.close()


class SplitDNA:
    def __init__(self, args):
        self.fq1 = args.fq1.split(",")
        self.fq2 = args.fq2.split(",")
        self.bcumi_runner = parse_chemistry.BcUmi("AccuraSCOPE_DNA", strict=True)
        _pattern_dict, bc_file = parse_chemistry.get_pattern_dict_and_bc("AccuraSCOPE_DNA")
        self.barcode_name = get_barcode_name(bc_file[0], args.well_name)
        self.outdir = Path("./dna")

    @utils.add_log
    def run(self):
        out_handles_r1 = {}
        out_handles_r2 = {}
        total_reads = 0
        valid_reads = 0
        signal_reads = 0
        name_read_count = defaultdict(int)
        namesheet_df_lines = []
        BC_LEN = 8

        # mkdir
        fastq_dir = self.outdir / "fastqs"
        if not fastq_dir.exists():
            fastq_dir.mkdir(parents=True, exist_ok=True)
        for fq1, fq2 in zip(self.fq1, self.fq2):
            with pysam.FastxFile(fq1) as handle1, pysam.FastxFile(fq2) as handle2:
                for read1, read2 in zip(handle1, handle2):
                    total_reads += 1
                    if total_reads % 1_000_000 == 0:
                        sys.stderr.write(f"Processed {total_reads} reads...\n")
                    valid, _corrected, corrected_seq, _umi = self.bcumi_runner.get_bc_umi(read1.sequence)
                    if not valid:
                        continue
                    valid_reads += 1
                    barcode = corrected_seq
                    name = self.barcode_name[barcode]
                    name_read_count[name] += 1

                    if not name.startswith("noise_"):
                        signal_reads += 1
                        if name not in out_handles_r1:
                            r1_file = fastq_dir / f"{name}_R1.fq.gz"
                            r2_file = fastq_dir / f"{name}_R2.fq.gz"
                            namesheet_df_lines.append(
                                {
                                    "patient": "patient1",
                                    "sample": name,
                                    "lane": 1,
                                    "fastq_1": str(r1_file.absolute()),
                                    "fastq_2": str(r2_file.absolute()),
                                }
                            )
                            out_handles_r1[name] = utils.generic_open(r1_file, "wt", compresslevel=1)
                            out_handles_r2[name] = utils.generic_open(r2_file, "wt", compresslevel=1)

                        out_handles_r1[name].write(utils.fastq_str(read1.name, read1.sequence[BC_LEN:], read1.quality))
                        out_handles_r2[name].write(utils.fastq_str(read2.name, read2.sequence, read2.quality))

        # output metrics
        metrics_file = self.outdir / "dna_metrics.txt"
        valid_reads_percent = round(valid_reads / total_reads * 100, 2) if total_reads > 0 else 0
        signal_reads_percent = round(signal_reads / total_reads * 100, 2) if total_reads > 0 else 0
        metrics_dict = {
            "total_reads": total_reads,
            "valid_reads": f"{valid_reads} ({valid_reads_percent}%)",
            "signal_reads": f"{signal_reads} ({signal_reads_percent}%)",
        }
        with open(metrics_file, "w") as mf:
            for k, v in metrics_dict.items():
                mf.write(f"{k}: {v}\n")

        # read count per name
        df = pd.DataFrame.from_dict(name_read_count, orient="index", columns=["read_count"])
        df.sort_values(by=["read_count"], ascending=False, inplace=True)
        name_count_file = self.outdir / "dna_read_count.txt"
        df.to_csv(name_count_file, index_label="name", sep="\t")

        # namesheet
        namesheet_file = self.outdir / "dna_samplesheet.csv"
        namesheet_df = pd.DataFrame(namesheet_df_lines)
        namesheet_df.to_csv(namesheet_file, index=False)

        # close file handles
        for handle in out_handles_r1.values():
            handle.close()
        for handle in out_handles_r2.values():
            handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="Split AccuraSCOPE RNA or DNA fastq into well fastq based on well name mapping",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--fq1",
        required=True,
        help="Input read1 fastq file (R1)",
    )
    parser.add_argument(
        "--fq2",
        required=True,
        help="Input read2 fastq file (R2)",
    )
    parser.add_argument(
        "--well_name",
        required=True,
        help="Well name file, tab delimited, two columns: well number and user defined name(e,g. 1\tcell_1), no header",
    )
    parser.add_argument(
        "--library_type",
        required=True,
        choices=["rna", "dna"],
        help="Library type, rna or dna",
    )
    args = parser.parse_args()
    if args.library_type == "rna":
        SplitRNA(args).run()
    elif args.library_type == "dna":
        SplitDNA(args).run()


if __name__ == "__main__":
    main()
