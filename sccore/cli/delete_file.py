#!/usr/bin/env python3
import os
import sys
import argparse
from multiprocessing import Pool, cpu_count


CURRENT_UID = os.getuid()


def delete_by_suffix_recursive(base_dir, suffixes, dry_run=False):
    try:
        with os.scandir(base_dir) as it:
            for entry in it:
                path = entry.path

                if entry.is_dir(follow_symlinks=False):
                    delete_by_suffix_recursive(path, suffixes, dry_run)

                elif entry.is_file(follow_symlinks=False):
                    if entry.name.endswith(suffixes):
                        if dry_run:
                            print(f"[DRY-RUN] {path}")
                        else:
                            try:
                                os.remove(path)
                                print(f"Deleted: {path}", flush=True)
                            except PermissionError:
                                pass
                            except FileNotFoundError:
                                pass

    except PermissionError:
        return
    except FileNotFoundError:
        return


def worker(args):
    path, suffixes, dry_run, check_uid = args

    if check_uid:
        try:
            if os.stat(path).st_uid != CURRENT_UID:
                # 顶层目录不是当前用户，直接跳过
                print(f"[SKIP UID] {path}")
                return
        except PermissionError:
            return
        except FileNotFoundError:
            return

    print(f"Processing: {path}")
    delete_by_suffix_recursive(path, suffixes, dry_run)


def parallel_delete(base_dir, suffixes, dry_run=False, workers=None, check_uid=False):
    tasks = []

    try:
        with os.scandir(base_dir) as it:
            for entry in it:
                path = entry.path

                if entry.is_dir(follow_symlinks=False):
                    tasks.append((path, suffixes, dry_run, check_uid))

    except PermissionError:
        print(f"[WARN] No permission: {base_dir}", file=sys.stderr)
        return

    if not tasks:
        return

    if workers is None:
        workers = min(len(tasks), cpu_count())

    print(f"Starting {workers} worker(s) under {base_dir} " f"(check_uid={check_uid})")

    with Pool(processes=workers) as pool:
        pool.map(worker, tasks)


SUFFIX_DEFAULT = [".bam", ".fq", ".fastq", ".bdg", "150bp.tsv"]


def main():
    parser = argparse.ArgumentParser(description="Parallel recursive delete files by suffix (permission-safe)")
    parser.add_argument("base_dir", help="Base directory to scan")
    parser.add_argument(
        "-s", "--suffix", default=",".join(SUFFIX_DEFAULT), help="Comma-separated suffixes, e.g. .bam,.fastq"
    )
    parser.add_argument(
        "-n", "--dry-run", action="store_true", help="Dry run: print files to be deleted without deleting"
    )
    parser.add_argument("-j", "--jobs", type=int, default=None, help="Number of worker processes")
    parser.add_argument(
        "--check_uid", action="store_true", help="Skip top-level subdirectories not owned by current user"
    )

    args = parser.parse_args()

    suffixes = tuple(s.strip() for s in args.suffix.split(",") if s.strip())

    if not suffixes:
        print("Error: no valid suffix provided", file=sys.stderr)
        sys.exit(1)

    parallel_delete(
        base_dir=args.base_dir,
        suffixes=suffixes,
        dry_run=args.dry_run,
        workers=args.jobs,
        check_uid=args.check_uid,
    )


if __name__ == "__main__":
    main()
