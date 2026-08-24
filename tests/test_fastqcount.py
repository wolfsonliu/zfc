import gzip
import re
from collections import Counter

import pytest

from zfc.fastqcount import (
    _count,
    _iter_batches,
    _iter_paired_records,
    _load_library,
    _match_barcode_pair,
    _match_guide,
    _open_text,
    _write_output,
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


def test_iter_paired_records_yields_sequences(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, ["AAA", "CCC"])
    _write_fastq(r2, ["TTT", "GGG"])
    with open(str(r1)) as h1, open(str(r2)) as h2:
        pairs = list(_iter_paired_records(h1, h2))
    assert pairs == [("AAA", "TTT"), ("CCC", "GGG")]


def test_iter_paired_records_mismatch_raises(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, ["AAA", "CCC"])
    _write_fastq(r2, ["TTT"])
    with open(str(r1)) as h1, open(str(r2)) as h2:
        with pytest.raises(ValueError, match="unequal"):
            list(_iter_paired_records(h1, h2))


def test_iter_paired_records_truncated_raises(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, ["AAA"])
    r2.write_text("@read1\nTTT\n+\n")  # missing quality line (truncated)
    with open(str(r1)) as h1, open(str(r2)) as h2:
        with pytest.raises(ValueError, match="truncated"):
            list(_iter_paired_records(h1, h2))


def test_iter_paired_records_empty(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    r1.write_text("")
    r2.write_text("")
    with open(str(r1)) as h1, open(str(r2)) as h2:
        assert list(_iter_paired_records(h1, h2)) == []


def test_iter_batches_splits_and_partial(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, ["AAA", "CCC", "GGG"])
    _write_fastq(r2, ["TTT", "AAA", "CCC"])
    batches = list(_iter_batches(str(r1), str(r2), 2))
    assert [len(batch) for batch in batches] == [2, 1]
    assert batches[0] == [("AAA", "TTT"), ("CCC", "AAA")]
    assert batches[1] == [("GGG", "CCC")]


def test_iter_batches_reads_gzip(tmp_path):
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    _write_gzip_fastq(r1, ["AAA", "CCC"])
    _write_gzip_fastq(r2, ["TTT", "GGG"])
    batches = list(_iter_batches(str(r1), str(r2), 2))
    assert batches == [[("AAA", "TTT"), ("CCC", "GGG")]]


def test_match_guide_basic():
    pattern = re.compile(GUIDE_PATTERN)
    assert _match_guide(_guide_seq(GUIDE_A), "T" * 20, pattern) == GUIDE_A


def test_match_guide_via_reverse_complement():
    pattern = re.compile(GUIDE_PATTERN)
    forward = _guide_seq(GUIDE_A)
    seq1 = reverse_complement(forward)
    assert _match_guide(seq1, "T" * 20, pattern) == GUIDE_A


def test_match_barcode_pair_anchor_on_read1():
    anchor = re.compile(ANCHOR_PATTERN)
    barcode_pattern = re.compile(BARCODE_PATTERN)
    seq1 = _anchor_seq(GUIDE_A, BARCODE_G)
    mate_struct = reverse_complement("TGGA" + BARCODE_G + "AACA")
    seq2 = "NN" + mate_struct + "NN"
    assert _match_barcode_pair(seq1, seq2, anchor, barcode_pattern) == (
        GUIDE_A, BARCODE_G, BARCODE_G
    )


def test_count_guide_raw_serial(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    seqs1 = [_guide_seq(GUIDE_A), _guide_seq(GUIDE_C), _guide_seq(GUIDE_A)]
    _write_fastq(r1, seqs1)
    _write_fastq(r2, ["T" * 20, "T" * 20, "T" * 20])
    counts, pairs = _count(
        "guide", {"guide": GUIDE_PATTERN}, None, str(r1), str(r2), 1, 10
    )
    assert pairs == 3
    assert counts == Counter({GUIDE_A: 2, GUIDE_C: 1})


def test_count_guide_with_library_counts_only_known_guides(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    lib = tmp_path / "lib.txt"
    _write_fastq(r1, [_guide_seq(GUIDE_A), _guide_seq(GUIDE_A)])
    _write_fastq(r2, ["T" * 20, "T" * 20])
    lib.write_text("gene\tguide\nGeneX\t{}\nGeneY\t{}\n".format(GUIDE_A, GUIDE_C))
    _entries, index = _load_library(str(lib), "guide")
    counts, pairs = _count(
        "guide", {"guide": GUIDE_PATTERN}, index, str(r1), str(r2), 1, 10
    )
    assert pairs == 2
    assert counts == Counter({("GeneX", GUIDE_A): 2})


def test_write_output_guide_raw(tmp_path):
    out = tmp_path / "out.tsv"
    _write_output(Counter({GUIDE_A: 2, GUIDE_C: 1}), None, "guide", str(out))
    assert out.read_text() == "guide\tcount\n" + GUIDE_A + "\t2\n" + GUIDE_C + "\t1\n"


def test_write_output_guide_zfc_ready_zero_fill(tmp_path):
    out = tmp_path / "out.tsv"
    entries = [("GeneX", GUIDE_A), ("GeneY", GUIDE_C)]
    counts = Counter({("GeneX", GUIDE_A): 5})
    _write_output(counts, entries, "guide", str(out))
    expected = (
        "gene\tguide\tbarcode\tcount\n"
        "GeneX\t{0}\t{0}\t5\n"
        "GeneY\t{1}\t{1}\t0\n"
    ).format(GUIDE_A, GUIDE_C)
    assert out.read_text() == expected


def test_write_output_barcode_raw(tmp_path):
    out = tmp_path / "out.tsv"
    _write_output(Counter({(GUIDE_A, BARCODE_G): 3}), None, "barcode", str(out))
    assert out.read_text() == (
        "guide\tbarcode\tcount\n{0}\t{1}\t3\n".format(GUIDE_A, BARCODE_G)
    )


def test_write_output_barcode_zfc_ready_zero_fill(tmp_path):
    out = tmp_path / "out.tsv"
    other = "T" * 6
    entries = [("GeneX", GUIDE_A, BARCODE_G), ("GeneY", GUIDE_A, other)]
    counts = Counter({("GeneX", GUIDE_A, BARCODE_G): 7})
    _write_output(counts, entries, "barcode", str(out))
    expected = (
        "gene\tguide\tbarcode\tcount\n"
        "GeneX\t{0}\t{1}\t7\n"
        "GeneY\t{0}\t{2}\t0\n"
    ).format(GUIDE_A, BARCODE_G, other)
    assert out.read_text() == expected


def test_write_output_stdout_branch(capsys):
    _write_output(Counter({GUIDE_A: 1}), None, "guide", "-")
    captured = capsys.readouterr()
    assert captured.out == "guide\tcount\n{0}\t1\n".format(GUIDE_A)