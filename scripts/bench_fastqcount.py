"""Benchmark zfc-count throughput on paired-end FASTQ."""

import argparse
import os
import time

from zfc.fastqcount import _count

DEFAULT_GUIDE = "ACCG([ATGC]{20})GTTT"
DEFAULT_ANCHOR = "ACCG([ATGC]{20})GTTT[ATGC]{5,35}TGGA([ATCG]{6})AACA"
DEFAULT_BARCODE = "TGGA([ATCG]{6})AACA"


def main(argv=None):
    parser = argparse.ArgumentParser(description="Benchmark zfc-count throughput.")
    parser.add_argument("-1", "--read1", required=True)
    parser.add_argument("-2", "--read2", required=True)
    parser.add_argument("--mode", choices=["guide", "barcode"], default="guide")
    parser.add_argument("--workers", type=int, default=0)
    args = parser.parse_args(argv)

    patterns = (
        {"guide": DEFAULT_GUIDE}
        if args.mode == "guide"
        else {"anchor": DEFAULT_ANCHOR, "barcode": DEFAULT_BARCODE}
    )
    workers = args.workers if args.workers > 0 else min(os.cpu_count() or 1, 8)

    start = time.perf_counter()
    counts, pairs = _count(
        args.mode, patterns, None, args.read1, args.read2, workers, 50000
    )
    elapsed = time.perf_counter() - start
    print(
        "mode={} workers={} pairs={} distinct_keys={} elapsed={:.2f}s rate={:.0f} pairs/s".format(
            args.mode, workers, pairs, len(counts), elapsed, pairs / elapsed
        )
    )


if __name__ == "__main__":
    main()