#!usr/bin/env python3
"""
create github test data
1. run SamplesheetGenerator to get sample_prefix_pair_fq and sample_barcode
2. run run() to get test data
3. generate tess data github samplesheet.csv
"""

import argparse
from collections import defaultdict
import subprocess
import os

from sccore.cli.manifest import SamplesheetGenerator
import pandas as pd

WEBSITE = "github"
ORG = "singleron-RD"


class GithubSamplesheet(SamplesheetGenerator):
    def __init__(self, manifest_file, folders, match, assay, dir_name):
        super().__init__(manifest_file, folders, match)
        self.prefix = f"https://{WEBSITE}.com/{ORG}/{assay}_test_data/raw/master/"

    def find_files(self):
        """
        Recursively search the specified folders for fastq files and (optionally) matched barcode files.
        """
        samples = set(self.manifest.values())
        for folder in self.folders:
            for root, _, files in os.walk(folder):
                print(root, files)
                for fn in files:
                    if not fn.endswith(".gz"):
                        continue
                    abspath = os.path.join(self.prefix, root, fn)
                    if self.match:
                        sample = self.get_barcode(root, fn, samples)
                        if sample:
                            self.sample_barcode[sample] = abspath
                            continue
                    sample, prefix, pair = self.get_fastq(fn, self.manifest)
                    if sample:
                        self.sample_prefix_pair_fq[sample][prefix][pair].append(abspath)
                        self.max_pair = max(self.max_pair, pair)
                        self.n_fq += 1
            print(f"Found {self.n_fq} fastqs in folder {folder}")
            if self.match:
                print(f"Found {len(self.sample_barcode)} match barcodes.tsv.gz in folder {folder}")


class CreateTestData:
    def __init__(self, pair_fq, match_barcode, dir_name, assay, n_reads, skip=1e5):
        self.pair_fq = pair_fq
        self.match_barcode = match_barcode
        self.n_reads = int(n_reads)
        self.skip = int(skip)
        self.dir_name = dir_name
        self.assay = assay

        self.manifest = {}
        self.manifest_fn = f"{dir_name}/manifest.csv"

    def make_folder(self):
        if not os.path.exists(self.dir_name):
            os.mkdir(self.dir_name)
        if self.match_barcode:
            for sample in ["X", "Y"]:
                matrix_dir = f"{self.dir_name}/match_barcode/{sample}.matrix/filtered"
                if not os.path.exists(matrix_dir):
                    os.makedirs(matrix_dir)
                # copy match_barcode to matrix_dir
                cmd = f"cp {self.match_barcode} {matrix_dir}"
                print(cmd)
                subprocess.check_call(cmd, shell=True)

    def subset_data(self, sample, repeat):
        for pair in self.pair_fq:
            fq = self.pair_fq[pair][0]
            for r in range(1, repeat + 1):
                total = (self.skip + self.n_reads * r) * 4
                prefix = f"prefix{sample}"
                fq_name = f"{prefix}_00{r}_R{pair}.fq.gz"
                self.manifest[prefix] = sample
                cmd = f"zcat {fq} | head -n {total} | tail -n {self.n_reads*4} | gzip > {self.dir_name}/{fq_name}"
                print(cmd)
                subprocess.check_call(cmd, shell=True)

    def write_manifest(self):
        dic = defaultdict(list)
        for prefix in self.manifest:
            sample = self.manifest[prefix]
            dic["sample"].append(sample)
            dic["prefix"].append(prefix)
        df = pd.DataFrame(dic)
        df.to_csv(self.manifest_fn, index=False)

    def write_samplesheet(self):
        match = True if self.match_barcode else False
        runner = GithubSamplesheet(self.manifest_fn, [self.dir_name], match, self.assay, self.dir_name)
        runner.run()
        samplesheet_fn = f"{self.dir_name}/samplesheet.csv"
        # move to dir_name
        cmd = f"mv samplesheet.csv {samplesheet_fn}"
        subprocess.check_call(cmd, shell=True)
        print(f"{runner.prefix}{samplesheet_fn}")

    def run(self):
        self.make_folder()
        self.subset_data("X", 2)
        self.subset_data("Y", 1)
        self.write_manifest()
        self.write_samplesheet()


def main():
    parser = argparse.ArgumentParser()
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
    parser.add_argument("--assay", "-a", help="assay name", required=True)
    parser.add_argument("--dir_name", "-d", help="dir name", required=True)
    parser.add_argument("--n_reads", "-n", help="number of reads in each fastq", default=5 * 1e4, type=int)
    # add dry_run
    parser.add_argument("--dry_run", help="dry run", action="store_true")
    args = parser.parse_args()

    # get sample_prefix_pair_fq and sample_barcode
    folders = args.folders.split(",")
    ssg = SamplesheetGenerator(
        args.manifest,
        folders,
        args.match,
    )
    ssg.find_files()
    sample = list(ssg.sample_prefix_pair_fq.keys())[0]
    prefix = list(ssg.sample_prefix_pair_fq[sample].keys())[0]
    pair_fq = ssg.sample_prefix_pair_fq[sample][prefix]
    match_barcode = ssg.sample_barcode.get(sample, None)

    runner = CreateTestData(pair_fq, match_barcode, args.dir_name, args.assay, args.n_reads)
    runner.run()


if __name__ == "__main__":
    main()
