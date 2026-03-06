import argparse
import random
import sys


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--length", "-l", help="length of barcode", type=int, required=True)
    parser.add_argument("--min-dist", "-d", help="minimum edit distance", type=int, required=True)
    parser.add_argument("--count", "-c", help="number of barcodes", type=int, required=True)
    parser.add_argument("--max-stretch", "-s", help="max stretch of the same nucleotide", type=int, required=False)
    parser.add_argument("--gc_min", help="gc min fraction", default=0.4)
    parser.add_argument("--gc_max", help="gc max fraction", default=0.6)
    parser.add_argument("--barcode_exist", help="exited barcode file", default=None)

    return parser.parse_args()


def gc_pass(s, gc_min, gc_max):
    n = len(s)
    gc_count = s.count("G") + s.count("C")
    if gc_count < n * float(gc_min) or gc_count > n * float(gc_max):
        return False
    return True


def hamming(a, b):
    """hamming distance between two strings"""
    dist = 0
    for i, j in zip(a, b):
        if i != j:
            dist += 1
    return dist


def max_stretch(s):
    """returns max stretch of the same letter in string"""
    max_stretch = 0
    last = None
    for n in s:
        if last is None:
            last = n
            stretch = 0
            continue
        if n == last:
            stretch += 1
            if stretch > max_stretch:
                max_stretch = stretch
        else:
            stretch = 0
        last = n
    return max_stretch


def random_nucleotide_sequence(length, alphabet="GATC"):
    """generate a random nucleotide sequence"""
    return "".join(random.sample(alphabet, 1)[0] for i in range(0, length))


def one_col_to_list(path) -> list[str]:
    with open(path) as f:
        return [line.strip() for line in f.readlines()]


def generate_barcodes(barcodes, **kwargs):
    """generate barcodes given constrains"""

    # spinner = cycle('|/-\\')

    index = 0
    failed_index = 0

    while True:
        if len(barcodes) == kwargs["count"]:
            break

        new_barcode = random_nucleotide_sequence(kwargs["length"])

        keep = True

        if new_barcode in barcodes:
            keep = False
            continue

        if keep:
            keep = gc_pass(new_barcode, kwargs["gc_min"], kwargs["gc_max"])

        if kwargs.get("max_stretch", False):
            if max_stretch(new_barcode) >= kwargs["max_stretch"]:
                keep = False

        if keep:
            for barcode in barcodes:
                if hamming(new_barcode, barcode) < kwargs["min_dist"]:
                    keep = False
                    break

        if keep:
            index += 1
            sys.stderr.write(str(index) + " done\n")
            # sys.stderr.write(spinner.next())
            barcodes.append(new_barcode)
        else:
            failed_index += 1
            if failed_index % 10000 == 0:
                sys.stderr.write(str(failed_index) + " failed\n")

    sys.stderr.write("finished")

    # make sure there are no duplicates!
    assert len(set(barcodes)) == len(barcodes)

    return barcodes


def main():
    """the guts"""

    args = parse_args()

    barcodes = []
    if args.barcode_exist:
        barcodes = one_col_to_list(args.barcode_exist)

    sys.stderr.write("generating barcodes...")
    barcodes = generate_barcodes(
        barcodes,
        min_dist=args.min_dist,
        count=args.count,
        length=args.length,
        max_stretch=args.max_stretch,
        gc_min=float(args.gc_min),
        gc_max=float(args.gc_max),
    )

    with open(f"{args.length}bp_bc.txt", "wt") as f:
        for bc in barcodes:
            f.write(bc + "\n")


if __name__ == "__main__":
    main()
