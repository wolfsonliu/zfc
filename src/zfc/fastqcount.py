"""Count sgRNAs (and iBAR barcodes) from paired-end FASTQ files.

Output-compatible Python replacement for the repo bash tools
``bin/library_count_sgrna`` and ``bin/library_count_sgrna_with_barcode``.
"""

from __future__ import annotations

import argparse
import gzip
import multiprocessing as mp
import os
import re
import sys
from collections import Counter

COMPLEMENT = str.maketrans("ATGC", "TACG")

DEFAULT_GUIDE = "ACCG([ATGC]{20})GTTT"
DEFAULT_ANCHOR = "ACCG([ATGC]{20})GTTT[ATGC]{5,35}TGGA([ATCG]{6})AACA"
DEFAULT_BARCODE = "TGGA([ATCG]{6})AACA"


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA string (A<->T, G<->C)."""
    return seq.translate(COMPLEMENT)[::-1]


def _open_text(path: str):
    """Open a file as UTF-8 text, transparently decompressing gzip input."""
    with open(path, "rb") as handle:
        magic = handle.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(path, "rt")
    return open(path, "r")


def _load_library(path: str, mode: str):
    """Parse the library file into an entry list and a lookup index.

    The first line is skipped when it looks like a header (first field is
    ``gene``). ``index`` maps a guide (guide mode) or a (guide, barcode)
    (barcode mode) to a gene; the first occurrence wins on duplicates.
    """
    if mode not in ("guide", "barcode"):
        raise ValueError("unknown mode: {!r}".format(mode))
    entries = []
    index = {}
    with _open_text(path) as handle:
        for lineno, line in enumerate(handle, 1):
            line = line.rstrip("\r\n")
            if lineno == 1 and line.split("\t", 1)[0] == "gene":
                continue
            if not line:
                continue
            cols = line.split("\t")
            if mode == "guide":
                if len(cols) < 2:
                    raise ValueError(
                        "library line {}: expected 'gene<tab>guide'".format(lineno)
                    )
                key = (cols[0], cols[1])
                index.setdefault(cols[1], cols[0])
            else:
                if len(cols) < 3:
                    raise ValueError(
                        "library line {}: expected 'gene<tab>guide<tab>barcode'".format(lineno)
                    )
                key = (cols[0], cols[1], cols[2])
                index.setdefault((cols[1], cols[2]), cols[0])
            entries.append(key)
    return entries, index


def _iter_paired_records(handle1, handle2):
    """Yield ``(seq1, seq2)`` for paired FASTQ files in lockstep."""
    while True:
        name1 = handle1.readline()
        name2 = handle2.readline()
        if not name1 and not name2:
            return
        if not name1 or not name2:
            raise ValueError("paired FASTQ files have unequal record counts")
        rec1 = [name1, handle1.readline(), handle1.readline(), handle1.readline()]
        rec2 = [name2, handle2.readline(), handle2.readline(), handle2.readline()]
        if any(line == "" for line in rec1 + rec2):
            raise ValueError("truncated FASTQ record")
        yield rec1[1].rstrip("\r\n"), rec2[1].rstrip("\r\n")


def _iter_batches(read1: str, read2: str, chunk_size: int):
    """Yield lists of at most ``chunk_size`` (seq1, seq2) tuples."""
    with _open_text(read1) as handle1, _open_text(read2) as handle2:
        batch = []
        for seq1, seq2 in _iter_paired_records(handle1, handle2):
            batch.append((seq1, seq2))
            if len(batch) >= chunk_size:
                yield batch
                batch = []
        if batch:
            yield batch