from sccore import parse_protocol, matrix, utils
import pandas as pd
from typing import Union


def create_pattern_args(pattern: str) -> str:
    """Create starsolo args relate to pattern"""
    pattern_dict = parse_protocol.parse_pattern(pattern)
    if len(pattern_dict["U"]) != 1:
        raise ValueError(
            f"Error: Wrong pattern:{pattern}. \n Solution: fix pattern so that UMI only have 1 position.\n"
        )
    ul = pattern_dict["U"][0].start
    ur = pattern_dict["U"][0].stop
    umi_len = ur - ul

    if len(pattern_dict["C"]) == 1:
        solo_type = "CB_UMI_Simple"
        start, stop = pattern_dict["C"][0].start, pattern_dict["C"][0].stop
        cb_start = start + 1
        cb_len = stop - start
        umi_start = ul + 1
        cb_str = f"--soloCBstart {cb_start} --soloCBlen {cb_len} --soloCBmatchWLtype 1MM "
        umi_str = f"--soloUMIstart {umi_start} --soloUMIlen {umi_len} "
    else:
        solo_type = "CB_UMI_Complex"
        cb_pos = " ".join([f"0_{x.start}_0_{x.stop-1}" for x in pattern_dict["C"]])
        umi_pos = f"0_{ul}_0_{ur-1}"
        cb_str = f"--soloCBposition {cb_pos} --soloCBmatchWLtype EditDist_2 "
        umi_str = f"--soloUMIposition {umi_pos} --soloUMIlen {umi_len} "

    return " ".join([f"--soloType {solo_type} ", cb_str, umi_str])


def create_v3_pattern_args() -> str:
    """GEXSCOPE-V3 has a 0-3 bases offset at the begining."""
    linker1 = "ACGATG"
    linker2 = "CATAGT"
    bc = "N" * 9
    linker_len = 6
    bc2_start = 9 + linker_len
    pattern_args = (
        "--soloType CB_UMI_Complex "
        f"--soloCBposition 2_0_2_8 2_{bc2_start}_2_{bc2_start+8} 3_1_3_9 "
        "--soloUMIposition 3_10_3_21 "
        f"--soloAdapterSequence {bc}{linker1}{bc}{linker2} "
        "--soloAdapterMismatchesNmax 1 "
        "--soloCBmatchWLtype EditDist_2 "
    )
    return pattern_args


def create_whitelist_args(whitelist_str) -> str:
    if not whitelist_str:
        res = "None"
    else:
        res = whitelist_str.strip()
        # nextflow copy remote file to current folder, so only keep file name.
        if res.startswith("http"):
            res = res.split("/")[-1]
        if res.endswith(".gz"):
            res = f"<(gzip -cdf {res})"
    res = f" --soloCBwhitelist {res} "
    return res


def create_solo_args(
    pattern_args: str,
    whitelist_args: str,
    sample: str,
    fq1: str,
    fq2: str,
    genomeDir: str,
    soloCellFilter: str,
    runThreadN: Union[str, int],
    clip3pAdapterSeq: str,
    outFilterMatchNmin: Union[str, int],
    soloFeatures: str,
    outSAMattributes: str,
    extra_starsolo_args: str,
) -> str:
    """Create all starsolo args"""
    read_command = "zcat" if fq1.strip().endswith(".gz") else "cat"
    cmd = (
        "STAR \\\n"
        f"{pattern_args} \\\n"
        f"{whitelist_args} \\\n"
        f"--outFileNamePrefix {sample}_ \\\n"
        f"--readFilesIn {fq2} {fq1} \\\n"
        f"--readFilesCommand {read_command} \\\n"
        f"--genomeDir {genomeDir} \\\n"
        f"--soloCellFilter {soloCellFilter} \\\n"
        f"--runThreadN {runThreadN} \\\n"
        f"--clip3pAdapterSeq {clip3pAdapterSeq} \\\n"
        f"--outFilterMatchNmin {outFilterMatchNmin} \\\n"
        f"--soloFeatures {soloFeatures} \\\n"
        f"--outSAMattributes {outSAMattributes} \\\n"
        f"{extra_starsolo_args} \\\n"
        "--outSAMtype BAM SortedByCoordinate \\\n"
        "--soloCellReadStats Standard \\\n"
        "--soloBarcodeReadLength 0 \\\n"
    )
    return cmd


class SoloSummary:
    def __init__(self, filtered_matrix, cellReadsStats, summary_txt):
        self.filtered_matrix = filtered_matrix
        self.cellReadsStats = cellReadsStats
        self.summary_txt = summary_txt

        self.cbs = matrix.CountMatrix.read_barcodes(self.filtered_matrix)

        # outputs
        self.origin_metrics = {}
        self.metrics = {}
        # umi_count: index is barcode, with two columns "UMI" and "mark"
        self.umi_count = pd.DataFrame()

    def parse_cellReadsStats(self):
        """update metrics and umi_count"""
        df = pd.read_csv(self.cellReadsStats, sep="\t", header=0, index_col=0)
        df = df.iloc[1:,]  # skip first line cb not pass whitelist
        umi_count = df.loc[:, ["nUMIunique", "countedU"]]  # keep dataframe format
        umi_count.rename(columns={"nUMIunique": "UMI"}, inplace=True)
        umi_count.sort_values(by="UMI", ascending=False, inplace=True)
        umi_count["mark"] = "UB"
        for cb in self.cbs:
            umi_count.loc[cb, "mark"] = "CB"
        self.umi_count = umi_count

        df = df.loc[
            :,
            [
                "cbMatch",
                "cbPerfect",
                "genomeU",
                "genomeM",
                "exonic",
                "intronic",
                "exonicAS",
                "intronicAS",
                "countedU",
            ],
        ]
        metrics = df.sum().to_dict()
        self.origin_metrics.update(metrics)

    def parse_summary(self):
        """update metrics"""
        data = utils.csv2dict(self.summary_txt)
        self.origin_metrics.update(data)
