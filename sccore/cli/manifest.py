#!/usr/bin/env python
"""
manifest -m manifest.csv -f folder1,folder2
create a samplesheet.csv file.
"""

import argparse
import os
import csv
import pathlib
import pandas as pd

from sccore import utils
from sccore.__init__ import VERSION


class SamplesheetGenerator:
    def __init__(self, manifest_file, folders, match=False, output_file="samplesheet.csv"):
        self.manifest_file = manifest_file
        self.match = match
        self.folders = folders
        self.output_file = output_file
        self.manifest = self.get_manifest()
        self.sample_prefix_pair_fq = utils.nested_defaultdict(dim=3, val_type=list)
        self.sample_barcode = {}
        self.max_pair = 1
        self.n_fq = 0

    def get_manifest(self):
        """
        Read the manifest CSV file and return a dictionary of {prefix: sample}.
        """
        manifest = {}
        with open(self.manifest_file, "r") as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                manifest[row["prefix"]] = row["sample"]
        return manifest

    def get_pair(self, file_name):
        """
        Determine the read pair (R1, R2, etc.) of a fastq file based on the file name.
        """
        for x in (1, 2, 3):
            for y in (f"_{x}_", f"_{x}.", f"_R{x}_", f"_R{x}."):
                if y in file_name:
                    return x

    def find_files(self):
        """
        Recursively search the specified folders for fastq files and (optionally) matched barcode files.
        """
        samples = set(self.manifest.values())
        print("Searching...")
        for folder in self.folders:
            folder = os.path.abspath(folder)
            for root, _, files in os.walk(folder):
                for fn in files:
                    if not fn.endswith(".gz"):
                        continue
                    abspath = os.path.join(root, fn)
                    if self.match:
                        sample = self.get_barcode(root, fn, samples)
                        if sample:
                            self.sample_barcode[sample] = abspath
                            print("Found match_barcode file:", abspath)
                            continue
                    sample, prefix, pair = self.get_fastq(fn, self.manifest)
                    if sample:
                        print("Found fastq file:", abspath)
                        self.sample_prefix_pair_fq[sample][prefix][pair].append(abspath)
                        self.max_pair = max(self.max_pair, pair)
                        self.n_fq += 1
            print(f"Found {self.n_fq} fastqs in folder {folder}")
            if self.match:
                print(f"Found {len(self.sample_barcode)} match barcodes.tsv.gz in folder {folder}")

    def get_barcode(self, root, fn, samples):
        """
        Return the sample name if the file is a matched barcode file.
        """

        if fn != "barcodes.tsv.gz":
            return None
        path = pathlib.PurePath(root)
        if path.name != "filtered":
            return None
        sample_dir = path.parent.name
        for sample in samples:
            if sample in sample_dir:
                return sample

    def get_fastq(self, fn, manifest):
        """
        Return the sample name, prefix, and read pair if the file is a fastq file.
        """
        if not (fn.endswith(".fastq.gz") or fn.endswith(".fq.gz")):
            return None, None, None
        for prefix, sample in manifest.items():
            if prefix + "_" in fn:
                return sample, prefix, self.get_pair(fn)
        return None, None, None

    def write_samplesheet(self):
        """
        Write the samplesheet CSV file.
        """
        rows = []
        for prefix, sample in self.manifest.items():
            pair_fqs = self.sample_prefix_pair_fq[sample][prefix]
            for x in range(1, self.max_pair + 1):
                pair_fqs[x].sort()

            for i in range(len(pair_fqs[1])):
                row = {"sample": sample}
                for x in range(1, self.max_pair + 1):
                    row[f"fastq_{x}"] = pair_fqs[x][i]
                if self.match:
                    row["match_barcode"] = self.sample_barcode[sample]
                rows.append(row)

        df = pd.DataFrame.from_dict(rows, orient="columns")
        df.to_csv(self.output_file, index=False)

    def run(self):
        print(f"Found {len(self.manifest)} samples in {self.manifest_file}")
        if self.match:
            print("Also searching for filtered/barcode.tsv.gz files")
        self.find_files()
        self.write_samplesheet()
        print("Done")


def main():
    parser = argparse.ArgumentParser(description="Generate samplesheet from fastq files.")
    parser.add_argument(
        "-m", "--manifest", required=True, help="Path to the manifest CSV file containing prefix-sample mapping."
    )
    parser.add_argument("-b", "--match", action="store_true", help="Whether to search for matched barcode.tsv.gz files")
    parser.add_argument(
        "-f",
        "--folders",
        required=True,
        help="Comma-separated paths to folders to search for fastq files and optional matched_barcode files.",
    )
    parser.add_argument("-o", "--output", help="Output samplesheet file.", default="samplesheet.csv")
    parser.add_argument("-v", "--version", action="version", version=f"{VERSION}")

    args = parser.parse_args()
    folders = args.folders.split(",")

    generator = SamplesheetGenerator(args.manifest, folders, match=args.match, output_file=args.output)
    generator.run()


if __name__ == "__main__":
    main()
