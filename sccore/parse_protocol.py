import itertools
import json
import os
import re
import sys
from importlib import resources
from collections import defaultdict

import pysam
from sccore import utils

logger = utils.get_logger(__name__)


def parse_pattern(pattern: str, allowed: str = "CLUNT") -> dict[str, list[slice]]:
    """Parse a pattern string into a dictionary of slices.

    >>> pattern_dict = parse_pattern("C8L16C8L16C8L1U12T18")
    >>> pattern_dict['C']
    [slice(0, 8, None), slice(24, 32, None), slice(48, 56, None)]
    >>> pattern_dict['L']
    [slice(8, 24, None), slice(32, 48, None), slice(56, 57, None)]
    """
    if not pattern:
        raise ValueError("Pattern cannot be an empty string")

    pattern_slices: dict[str, list[slice]] = {}
    p = re.compile(r"([A-Z])(\d+)")  # Compile the regex

    start = 0
    for char, length in p.findall(pattern):
        if char not in allowed:
            raise ValueError(f"Invalid character '{char}' in pattern: {pattern}")
        if char not in pattern_slices:
            pattern_slices[char] = []
        end = start + int(length)
        pattern_slices[char].append(slice(start, end))
        start = end
    return pattern_slices


def create_mismatch_seqs(seq: str, max_mismatch=1, allowed_bases="ACGTN") -> set[str]:
    """Create all sequences within a specified number of mismatches from the input sequence.

    >>> answer = set(["TCG", "AAG", "ACC", "ATG", "ACT", "ACN", "GCG", "ANG", "ACA", "ACG", "CCG", "AGG", "NCG"])
    >>> seq_set = create_mismatch_seqs("ACG")
    >>> seq_set == answer
    True
    """
    if max_mismatch < 0:
        raise ValueError("max_mismatch must be non-negative")
    if max_mismatch > len(seq):
        raise ValueError(f"max_mismatch ({max_mismatch}) cannot be greater than the sequence length ({len(seq)})")

    result = set()
    for locs in itertools.combinations(range(len(seq)), max_mismatch):
        seq_locs = [list(allowed_bases) if i in locs else [base] for i, base in enumerate(seq)]
        result.update("".join(p) for p in itertools.product(*seq_locs))
    return result


def create_mismatch_origin_dict(origin_seqs: list, n_mismatch=1) -> dict[str, str]:
    """Create a dictionary mapping sequences with mismatches to their original sequences(in whitelist).

    >>> origin_seqs = ["AACGTGAT", "AAACATCG"]
    >>> mismatch_dict = create_mismatch_origin_dict(origin_seqs)
    >>> mismatch_dict["AACGTGAA"] == "AACGTGAT"
    True
    """
    result = {}
    for origin_seq in origin_seqs:
        origin_seq = origin_seq.strip()
        if origin_seq == "":
            continue
        for mismatch_seq in create_mismatch_seqs(origin_seq, n_mismatch):
            result[mismatch_seq] = origin_seq
    return result


def create_mismatch_origin_dicts_from_whitelists(whitelists: list, n_mismatch: int = 1) -> tuple[list, list]:
    """Returns origin whitelist list and mismatch dict list.

    >>> whitelists = [resources.files("sccore.protocols").joinpath("whitelist/GEXSCOPE-V1/bc.txt")]
    >>> raw_list, mismatch_list = create_mismatch_origin_dicts_from_whitelists(whitelists)
    >>> len(raw_list) == len(mismatch_list)
    True
    """
    raw_list, mismatch_list = [], []
    for f in whitelists:
        barcodes = utils.read_one_col(f)
        raw_list.append(set(barcodes))
        barcode_mismatch_dict = create_mismatch_origin_dict(barcodes, n_mismatch)
        mismatch_list.append(barcode_mismatch_dict)

    return raw_list, mismatch_list


def check_seq_mismatch(seq_list, raw_list, mismatch_list):
    """
    Returns
        valid: True if seq in mismatch_list
        corrected: True if seq in mismatch_list but not in raw_list
        res: joined seq

    >>> seq_list = ['ATA', 'AAT', 'ATA']
    >>> correct_set_list = [{'AAA'},{'AAA'},{'AAA'}]
    >>> mismatch_dict_list = [create_mismatch_origin_dict(['AAA'])] * 3

    >>> check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
    (True, True, 'AAA_AAA_AAA')

    >>> seq_list = ['AAA', 'AAA', 'AAA']
    >>> check_seq_mismatch(seq_list, correct_set_list, mismatch_dict_list)
    (True, False, 'AAA_AAA_AAA')
    """
    valid = True
    corrected = False
    res = []
    for index, seq in enumerate(seq_list):
        if seq not in raw_list[index]:
            if seq not in mismatch_list[index]:
                valid = False
                res.append("")
            else:
                corrected = True
                res.append(mismatch_list[index][seq])
        else:
            res.append(seq)

    return valid, corrected, "_".join(res)


def get_protocol_dict():
    """
    Return:
    protocol_dict. Key: protocol name, value: protocol dict

    >>> protocol_dict = get_protocol_dict()
    >>> protocol_dict["GEXSCOPE-MicroBead"]["pattern_dict"]
    {'C': [slice(0, 12, None)], 'U': [slice(12, 20, None)]}
    """
    json_file = str(resources.files("sccore.protocols").joinpath("protocols.json"))
    protocol_dict = json.load(open(json_file))
    whitelist_dir = str(resources.files("sccore.protocols").joinpath("whitelist"))
    # add folder prefix
    for protocol in protocol_dict:
        cur = protocol_dict[protocol]
        bc = cur.get("bc", [])
        linker = cur.get("linker", [])
        if bc:
            cur["bc"] = [os.path.join(whitelist_dir, protocol, x) for x in bc]
        if linker:
            cur["linker"] = [os.path.join(whitelist_dir, protocol, x) for x in linker]
        cur["pattern_dict"] = parse_pattern(cur["pattern"])
    return protocol_dict


class Auto:
    """
    Auto detect singleron protocols from R1-read
    """

    def __init__(self, fq1_list, protocol_dict, max_read=10000):
        """
        Returns:
            protocol, protocol_dict[protocol]
        """
        self.fq1_list = fq1_list
        self.max_read = max_read
        self.protocol_dict = protocol_dict
        self.mismatch_dict = {}
        for protocol in self.protocol_dict:
            if "bc" in self.protocol_dict[protocol]:
                self.mismatch_dict[protocol] = create_mismatch_origin_dicts_from_whitelists(
                    self.protocol_dict[protocol]["bc"], 1
                )

    def run(self):
        """
        Returns:
            protocol, protocol_dict[protocol]
        """
        protocol = self.get_protocol()
        return protocol, self.protocol_dict[protocol]

    def get_protocol(self):
        """check protocol in the fq1_list"""
        fq_protocol = {}
        for fastq1 in self.fq1_list:
            protocol = self.get_fq_protocol(fastq1)
            fq_protocol[fastq1] = protocol
        if len(set(fq_protocol.values())) != 1:
            sys.exit(
                f"Error: multiple protocols are not allowed for one sample: {self.fq1_list}! \n" + str(fq_protocol)
            )
        protocol = list(fq_protocol.values())[0]
        return protocol

    def is_protocol(self, seq, protocol):
        """check if seq matches the barcode of protocol"""
        raw_list, mismatch_list = self.mismatch_dict[protocol]
        bc_list = [seq[x] for x in self.protocol_dict[protocol]["pattern_dict"]["C"]]
        valid, _corrected, _res = check_seq_mismatch(bc_list, raw_list, mismatch_list)
        return valid

    def seq_protocol(self, seq):
        """check if seq matches any protocol in protocol_dict"""
        for protocol in self.protocol_dict:
            if self.is_protocol(seq, protocol):
                return protocol
        return None

    def get_fq_protocol(self, fq1):
        protocol_readcount = defaultdict(int)

        fq = pysam.FastxFile(fq1)
        n = 0
        for read in fq:
            seq = read.sequence
            n += 1
            protocol = self.seq_protocol(seq)
            if protocol:
                protocol_readcount[protocol] += 1
            if n == self.max_read:
                break
        sorted_counts = sorted(protocol_readcount.items(), key=lambda x: x[1], reverse=True)
        logger.info(sorted_counts)

        protocol, read_counts = sorted_counts[0]
        percent = float(read_counts) / n
        if percent < 0.1:
            logger.error("Valid protocol read counts percent < 0.1")
            raise Exception("Auto protocol detection failed! ")
        elif percent < 0.5:
            logger.warning("Valid protocol read counts percent < 0.5")
        logger.info(f"{fq1}: {protocol}")

        return protocol


class AutoRNA(Auto):
    def __init__(self, fq1_list, max_read=10000):
        super().__init__(fq1_list, get_protocol_dict(), max_read)
        self.v3_linker_mismatch = create_mismatch_origin_dicts_from_whitelists(
            self.protocol_dict["GEXSCOPE-V3"]["linker"], 1
        )

    def v3_offset(self, seq):
        """
        return -1 if not v3

        >>> seq = "AT" + "TCGACTGTC" + "ACGATG" + "TTCTAGGAT" + "CATAGT" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner = AutoRNA([], "fake_sample")
        >>> runner.v3_offset(seq)
        2
        >>> seq = "TCGACTGTC" + "ACGATG" + "TTCTAGGAT" + "CATAGT" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.v3_offset(seq)
        0
        >>> seq = "TCGACTGTC" + "ATATAT" + "TTCTAGGAT" + "CATAGT" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.v3_offset(seq)
        -1
        """
        bc_len = 9
        linker_len = 6
        max_offset_len = 3 + 1  # allow for extra 1 bases
        for offset in range(max_offset_len + 1):
            first_linker_start = offset + bc_len
            second_linker_start = first_linker_start + linker_len + bc_len
            first_linker_seq = seq[first_linker_start : first_linker_start + linker_len]
            second_linker_seq = seq[second_linker_start : second_linker_start + linker_len]
            valid, _, _ = check_seq_mismatch([first_linker_seq, second_linker_seq], *self.v3_linker_mismatch)
            if valid:
                return offset
        return -1

    def seq_protocol(self, seq):
        """
        Returns: protocol or None

        >>> import tempfile
        >>> runner = AutoRNA([], "fake_sample")
        >>> seq = "AT" + "TCGACTGTC" + "ACGATG" + "TTCTAGGAT" + "CATAGT" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V3'
        >>> seq = "TCGACTGTC" + "ATCCACGTGCTTGAGA" + "TTCTAGGAT" + "TCAGCATGCGGCTACG" + "TGCACGAGA" + "C" + "CATATCAATGGG" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V2'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA" + "C" + "TCCGAAGCCCAT" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC" + "CTGTCT"
        >>> runner.seq_protocol(seq)
        'flv_rna'
        >>> seq = "NCAGATTC" + "TCGGTGACAGCCATAT" + "GTACGCAA" + "CGTAGTCAGAAGCTGA" + "CTGAGCCA"  + "TCCGAAGCC"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-V1'
        >>> seq = "ATCGATCGATCG" + "ATCGATCG" + "C" + "TTTTTTTTTT"
        >>> runner.seq_protocol(seq)
        'GEXSCOPE-MicroBead'
        """
        if self.v3_offset(seq) != -1:
            return "GEXSCOPE-V3"

        for protocol in ["GEXSCOPE-V2", "GEXSCOPE-V1", "AccuraCode"]:
            if self.is_protocol(seq, protocol):
                if protocol == "GEXSCOPE-V1":
                    if seq[self.protocol_dict["flv_rna"]["pattern_dict"]["L"][2]] == "CTGTCT":
                        return "flv_rna"
                return protocol

        # check if it is MicroBead
        if seq[16:20] != "TTTT" and seq[22:26] == "TTTT":
            return "GEXSCOPE-MicroBead"
