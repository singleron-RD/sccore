#!usr/bin/env python3
"""
Find barcode.tsv.gz under folders. Add underscore to each line and trim '-1' suffix.
"""

import argparse
import os

from sccore import utils

logger = utils.get_logger(__name__)


def format(bc):
    bc = bc.replace("-1", "")
    if len(bc) in (24, 27):
        sub = len(bc) // 3
        return "_".join([bc[:sub], bc[sub : 2 * sub], bc[2 * sub :]])


def need_format(bc):
    if "-1" in bc:
        return True
    if len(bc) in (24, 27) and "_" not in bc:
        return True
    return False


def run(folders):
    for folder in folders:
        for root, _, files in os.walk(folder):
            for file in files:
                if file == "barcodes.tsv.gz":
                    file_path = os.path.join(root, file)
                    bcs = utils.read_one_col(file_path)
                    if need_format(bcs[0]):
                        logger.info(f"Formatting {file_path}")
                    else:
                        logger.info(f"No need to format {file_path}")
                        continue
                    os.rename(file_path, file_path + ".origin")
                    bcs = [format(bc) for bc in bcs]
                    utils.write_one_col(bcs, file_path)


def main():
    parser = argparse.ArgumentParser(description="Add underscore to barcode.tsv.gz. Trim '-1' suffix.")
    parser.add_argument(
        "-f",
        "--folders",
        required=True,
        help="Comma-separated paths to folders to search for barcodes.tsv.gz",
    )

    args = parser.parse_args()
    folders = args.folders.split(",")
    run(folders)


if __name__ == "__main__":
    main()
