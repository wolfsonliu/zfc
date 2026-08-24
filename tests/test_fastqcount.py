import gzip
import re
from collections import Counter

import pytest

from zfc.fastqcount import (
    _open_text,
    reverse_complement,
)

GUIDE_PATTERN = "ACCG([ATGC]{20})GTTT"
ANCHOR_PATTERN = "ACCG([ATGC]{20})GTTT[ATGC]{5,35}TGGA([ATCG]{6})AACA"
BARCODE_PATTERN = "TGGA([ATCG]{6})AACA"

GUIDE_A = "A" * 20
GUIDE_C = "C" * 20
BARCODE_G = "G" * 6


def _write_fastq(path, sequences):
    lines = []
    for index, seq in enumerate(sequences, 1):
        lines.append("@read{}".format(index))
        lines.append(seq)
        lines.append("+")
        lines.append("I" * len(seq))
    path.write_text("\n".join(lines) + "\n")


def _write_gzip_fastq(path, sequences):
    lines = []
    for index, seq in enumerate(sequences, 1):
        lines.append("@read{}".format(index))
        lines.append(seq)
        lines.append("+")
        lines.append("I" * len(seq))
    with gzip.open(path, "wt") as handle:
        handle.write("\n".join(lines) + "\n")


def _guide_seq(guide):
    return "NNACCG" + guide + "GTTTNN"


def _anchor_seq(guide, barcode):
    return "NNACCG" + guide + "GTTTAAAAATGGA" + barcode + "AACANN"


def test_reverse_complement():
    assert reverse_complement("ATGC") == "GCAT"
    assert reverse_complement("AACCGGTT") == "AACCGGTT"


def test_open_text_detects_gzip(tmp_path):
    plain = tmp_path / "plain.fastq"
    plain.write_text("ACGT\n")
    gzed = tmp_path / "reads.fastq.gz"
    _write_gzip_fastq(gzed, ["ACGT"])
    with _open_text(str(plain)) as handle:
        assert handle.read() == "ACGT\n"
    with _open_text(str(gzed)) as handle:
        assert handle.read() == "@read1\nACGT\n+\nIIII\n"