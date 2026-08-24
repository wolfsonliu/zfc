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