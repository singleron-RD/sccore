#!usr/bin/env python3

import argparse
import pyfastx
import pysam
from sccore import parse_protocol, utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fq1")
    parser.add_argument("--bam")
    parser.add_argument("--assets_dir")
    parser.add_argument("--num", default=10**4, type=int)
    args = parser.parse_args()
    if not (args.fq1 or args.bam):
        raise ValueError("Please provide fq1 or bam")

    if args.bam:
        cnt = 0
        read_names = set()
        fh = pysam.AlignmentFile(args.bam, "rb")
        for read in fh:
            bc = read.get_tag("CB")
            if bc == "-" and read.query_name not in read_names:
                cnt += 1
                read_names.add(read.query_name)
                if cnt == args.num:
                    break
        print(f"{cnt} Read names extracted")

        fq1 = pyfastx.Fastx(args.fq1)
        found = 0
        invalid_fastq = open("invalid.txt", "wt")
        for name1, seq1, qual1 in fq1:
            if name1 in read_names:
                found += 1
                invalid_fastq.write(seq1 + "\n")
                if found == args.num:
                    break

    elif args.fq1:
        protocol_dict = parse_protocol.get_protocol_dict(args.assets_dir)
        v2_dict = protocol_dict["GEXSCOPE-V2"]
        v2_raw, v2_mismatch = parse_protocol.create_mismatch_origin_dicts_from_whitelists(v2_dict["bc"], 1)

        invalid_fastq = open("invalid.fastq", "wt")
        fq1 = pyfastx.Fastx(args.fq1)
        raw = valid_reads = invalid_reads = 0
        for name1, seq1, qual1 in fq1:
            raw += 1
            bc_list = [seq1[x] for x in v2_dict["pattern_dict"]["C"]]
            valid, _corrected, res = parse_protocol.check_seq_mismatch(bc_list, v2_raw, v2_mismatch)
            if not valid:
                invalid_fastq.write(utils.fastq_str(name1, seq1, qual1))
                invalid_reads += 1
            else:
                valid_reads += 1

            if raw == args.num:
                break
        print(f"Total reads: {raw}")
        print(f"Valid reads: {valid_reads}")
        print(f"Invalid reads: {invalid_reads}")


if __name__ == "__main__":
    main()
