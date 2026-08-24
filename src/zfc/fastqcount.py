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
    if chunk_size < 1:
        raise ValueError("chunk_size must be >= 1")
    with _open_text(read1) as handle1, _open_text(read2) as handle2:
        batch = []
        for seq1, seq2 in _iter_paired_records(handle1, handle2):
            batch.append((seq1, seq2))
            if len(batch) >= chunk_size:
                yield batch
                batch = []
        if batch:
            yield batch


def _match_guide(seq1, seq2, pattern):
    """Return the first guide captured across R1, R2, RC(R1), RC(R2)."""
    for seq in (seq1, seq2, reverse_complement(seq1), reverse_complement(seq2)):
        match = pattern.search(seq)
        if match:
            return match.group(1)
    return None


def _match_barcode_pair(seq1, seq2, anchor, barcode_pattern):
    """Return ``(guide, barcode1, barcode2)`` or None for a read pair.

    ``barcode1`` is captured by the anchor on the same read as the guide;
    ``barcode2`` is captured from the mate read by the barcode pattern.
    Candidates are tried in order: (R1, RC(R2)), (RC(R1), R2),
    (R2, RC(R1)), (RC(R2), R1).
    """
    candidates = (
        (seq1, reverse_complement(seq2)),
        (reverse_complement(seq1), seq2),
        (seq2, reverse_complement(seq1)),
        (reverse_complement(seq2), seq1),
    )
    for anchor_read, mate_read in candidates:
        match = anchor.search(anchor_read)
        if match:
            guide = match.group(1)
            barcode1 = match.group(2) if (match.lastindex or 0) >= 2 else None
            mate_match = barcode_pattern.search(mate_read)
            barcode2 = mate_match.group(1) if mate_match else None
            return guide, barcode1, barcode2
    return None


_WORKER_MODE = None
_WORKER_PATTERNS = None
_WORKER_LIBRARY = None


def _init_worker(mode, pattern_strings, library):
    global _WORKER_MODE, _WORKER_PATTERNS, _WORKER_LIBRARY
    _WORKER_MODE = mode
    _WORKER_PATTERNS = {name: re.compile(value) for name, value in pattern_strings.items()}
    _WORKER_LIBRARY = library


def _process_batch(batch):
    """Count one batch of read pairs. Returns (partial Counter, pair count)."""
    counts = Counter()
    mode = _WORKER_MODE
    library = _WORKER_LIBRARY
    patterns = _WORKER_PATTERNS
    if mode == "guide":
        pattern = patterns["guide"]
        for seq1, seq2 in batch:
            guide = _match_guide(seq1, seq2, pattern)
            if guide is None:
                continue
            if library is not None:
                if guide in library:
                    counts[(library[guide], guide)] += 1
            else:
                counts[guide] += 1
    else:
        anchor = patterns["anchor"]
        barcode_pattern = patterns["barcode"]
        for seq1, seq2 in batch:
            hit = _match_barcode_pair(seq1, seq2, anchor, barcode_pattern)
            if hit is None:
                continue
            guide, barcode1, barcode2 = hit
            if library is not None:
                if barcode1 is not None and (guide, barcode1) in library:
                    counts[(library[(guide, barcode1)], guide, barcode1)] += 1
                elif barcode2 is not None and (guide, barcode2) in library:
                    counts[(library[(guide, barcode2)], guide, barcode2)] += 1
            else:
                barcode = barcode1 if barcode1 is not None else barcode2
                counts[(guide, barcode)] += 1
    return counts, len(batch)


def _count(mode, pattern_strings, library, read1, read2, workers, chunk_size):
    """Merge partial counts from serial or pooled batch processing."""
    total = Counter()
    pairs = 0
    if workers <= 1:
        _init_worker(mode, pattern_strings, library)
        for batch in _iter_batches(read1, read2, chunk_size):
            partial, size = _process_batch(batch)
            total.update(partial)
            pairs += size
    else:
        with mp.Pool(
            workers,
            initializer=_init_worker,
            initargs=(mode, pattern_strings, library),
        ) as pool:
            for partial, size in pool.imap(
                _process_batch,
                _iter_batches(read1, read2, chunk_size),
                chunksize=1,
            ):
                total.update(partial)
                pairs += size
    return total, pairs


def _write_output(counts, entries, mode, output):
    """Write the headered count table to ``output`` ('-' means stdout)."""
    close = False
    if output in (None, "-"):
        handle = sys.stdout
    else:
        handle = open(output, "w")
        close = True
    try:
        if mode == "guide":
            if entries is not None:
                handle.write("gene\tguide\tbarcode\tcount\n")
                for gene, guide in sorted(set(entries)):
                    count = counts.get((gene, guide), 0)
                    handle.write("{}\t{}\t{}\t{}\n".format(gene, guide, guide, count))
            else:
                handle.write("guide\tcount\n")
                for guide in sorted(counts):
                    handle.write("{}\t{}\n".format(guide, counts[guide]))
        else:
            if entries is not None:
                handle.write("gene\tguide\tbarcode\tcount\n")
                for gene, guide, barcode in sorted(set(entries)):
                    count = counts.get((gene, guide, barcode), 0)
                    handle.write(
                        "{}\t{}\t{}\t{}\n".format(gene, guide, barcode, count)
                    )
            else:
                handle.write("guide\tbarcode\tcount\n")
                for guide, barcode in sorted(counts):
                    handle.write(
                        "{}\t{}\t{}\n".format(guide, barcode, counts[(guide, barcode)])
                    )
    finally:
        if close:
            handle.close()