import gzip
import re
from collections import Counter

import pytest

from zfc.fastqcount import (
    _load_library,
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


def test_load_library_guide_skips_header(tmp_path):
    path = tmp_path / "lib.txt"
    path.write_text("gene\tguide\nGeneX\t{}\nGeneY\t{}\n".format(GUIDE_A, GUIDE_C))
    entries, index = _load_library(str(path), "guide")
    assert entries == [("GeneX", GUIDE_A), ("GeneY", GUIDE_C)]
    assert index == {GUIDE_A: "GeneX", GUIDE_C: "GeneY"}


def test_load_library_barcode_skips_header(tmp_path):
    path = tmp_path / "lib.txt"
    path.write_text("gene\tguide\tbarcode\nGeneX\t{0}\t{1}\n".format(GUIDE_A, BARCODE_G))
    entries, index = _load_library(str(path), "barcode")
    assert entries == [("GeneX", GUIDE_A, BARCODE_G)]
    assert index == {(GUIDE_A, BARCODE_G): "GeneX"}


def test_load_library_no_header_first_line_is_data(tmp_path):
    path = tmp_path / "lib.txt"
    path.write_text("GeneX\t{}\n".format(GUIDE_A))
    entries, index = _load_library(str(path), "guide")
    assert entries == [("GeneX", GUIDE_A)]


def test_load_library_duplicate_guide_first_wins(tmp_path):
    path = tmp_path / "lib.txt"
    path.write_text("gene\tguide\nGeneX\t{0}\nGeneY\t{0}\n".format(GUIDE_A))
    entries, index = _load_library(str(path), "guide")
    assert entries == [("GeneX", GUIDE_A), ("GeneY", GUIDE_A)]
    assert index == {GUIDE_A: "GeneX"}


def test_load_library_reads_gzip(tmp_path):
    path = tmp_path / "lib.txt.gz"
    with gzip.open(path, "wt") as handle:
        handle.write("gene\tguide\nGeneX\t{}\n".format(GUIDE_A))
    entries, index = _load_library(str(path), "guide")
    assert entries == [("GeneX", GUIDE_A)]
    assert index == {GUIDE_A: "GeneX"}