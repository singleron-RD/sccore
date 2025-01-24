"""
extract barcode and UMI from fq1 and add them to the read name of fq2
"""

import pysam
from sccore import parse_protocol, utils


class Extract:
    def __init__(self, fq1_list, fq2_list, sample: str, protocol_dict, protocol):
        self.fq1_list = fq1_list
        self.fq2_list = fq2_list
        self.sample = sample
        self.protocol = protocol
        protocol_meta = protocol_dict[protocol]
        self.pattern_dict = protocol_meta["pattern_dict"]
        self.raw_list, self.mismatch_list = parse_protocol.create_mismatch_origin_dicts_from_whitelists(
            protocol_meta["bc"], 1
        )
        # v3
        if protocol == "GEXSCOPE-V3":
            self.offset_runner = parse_protocol.AutoRNA(fq1_list)
        # output
        self.out_fq = f"{sample}_R2.fq.gz"

    def get_bc_umi(self, seq):
        """
        >>> runner = Extract([],[],"fake", parse_protocol.get_protocol_dict(), "GEXSCOPE-V3")
        >>> seq = "AT" + "TCGACTGTC" + "ACGATG" + "TTCGAGGAT" + "CATAGT" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.get_bc_umi(seq)
        (True, False, 'TCGACTGTC_TTCGAGGAT_TGCACGAGA', 'CATATCAATGGG')
        """
        if self.protocol == "GEXSCOPE-V3":
            offset = self.offset_runner.v3_offset(seq)
            seq = seq[offset:]
        bc_list = [seq[x] for x in self.pattern_dict["C"]]
        valid, corrected, corrected_seq = parse_protocol.check_seq_mismatch(bc_list, self.raw_list, self.mismatch_list)
        if not valid:
            umi = None
        else:
            umi = seq[self.pattern_dict["U"][0]]
        return valid, corrected, corrected_seq, umi

    def run(self):
        raw_reads = valid_reads = corrected_reads = 0
        with utils.openfile(self.out_fq, "wt") as out_fh:
            for fq1, fq2 in zip(self.fq1_list, self.fq2_list):
                fq1 = pysam.FastxFile(fq1)
                fq2 = pysam.FastxFile(fq2)
                for e1, e2 in zip(fq1, fq2):
                    raw_reads += 1
                    valid, corrected, corrected_seq, umi = self.get_bc_umi(e1.sequence)
                    if valid:
                        valid_reads += 1
                        if corrected:
                            corrected_reads += 1
                        read_name = f"{corrected_seq}:{umi}:{raw_reads}"
                        out_fh.write(f"@{read_name}\n{e2.sequence}\n+\n{e2.quality}\n")  # type: ignore
        return raw_reads, valid_reads, corrected_reads
