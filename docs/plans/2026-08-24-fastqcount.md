# FASTQ 计数脚本 Python 重写（zfc-count）实现计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 用 Python 单遍流式重写两个 bash 计数脚本，产出 `zfc-count` CLI（支持 gzip、并行、`-l` 可选），输出满足 `zfc` 输入契约。

**Architecture:** 新增 `src/zfc/fastqcount.py`。主进程按 4 行记录边界流式读两个 FASTQ 并切块；worker（`ProcessPoolExecutor`）用预编译正则提取+计数，返回部分 `Counter`；主进程合并后用 broadcast 的 library 做基因映射并输出。serial（`--workers 1`）与 parallel 共用同一批处理函数。

**Tech Stack:** Python ≥3.11 标准库（`gzip`/`re`/`argparse`/`collections`/`multiprocessing`），pytest，ruff。**无新增运行时依赖。**

**设计文档:** `docs/specs/2026-08-24-fastqcount-design.md`

---

## 文件结构

- **创建** `src/zfc/fastqcount.py` — 计数核心：解析、匹配、计数、库映射、`main()`。
- **创建** `tests/test_fastqcount.py` — 单元 + CLI 端到端测试。
- **创建** `scripts/bench_fastqcount.py` — 吞吐基准脚本。
- **修改** `pyproject.toml` — 注册 `zfc-count = "zfc.fastqcount:main"`。
- **修改** `docs/data-format.md`、`docs/architecture.md`、（可选）`docs/development.md`、`README.md` — 标注新 CLI 与 bash 脚本 deprecated。

模块函数接口（本计划内固定，勿改名）：

| 函数 | 签名 | 职责 |
|---|---|---|
| `reverse_complement` | `(seq: str) -> str` | 反向互补 |
| `_open_text` | `(path: str) -> IO` | 按魔数返回 gz/plain 文本句柄 |
| `_require_groups` | `(pattern, name, minimum) -> re.Pattern` | 编译并校验捕获组数 |
| `_load_library` | `(path, mode) -> (entries, index)` | 解析 library（跳过 header） |
| `_iter_paired_records` | `(h1, h2) -> Iterator[(str, str)]` | 锁步产出 (seq1, seq2) |
| `_iter_batches` | `(r1, r2, chunk_size) -> Iterator[list]` | 分批 |
| `_match_guide` | `(seq1, seq2, pattern) -> str|None` | guide 模式提取 |
| `_match_barcode_pair` | `(seq1, seq2, anchor, barcode) -> tuple|None` | barcode 模式提取 |
| `_init_worker` | `(mode, pattern_strings, library) -> None` | worker 初始化 |
| `_process_batch` | `(batch) -> (Counter, int)` | 单批计数 |
| `_count` | `(...) -> (Counter, int)` | 编排串行/并行 |
| `_write_output` | `(counts, entries, mode, output) -> None` | 写表 |
| `_build_patterns` | `(args) -> dict` | 校验并组装 pattern dict |
| `main` | `(argv=None) -> None` | CLI |

---

## Task 1: 模块骨架 — `reverse_complement` 与 `_open_text`

**Files:**
- Create: `src/zfc/fastqcount.py`
- Create: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

创建 `tests/test_fastqcount.py`：

```python
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
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'zfc.fastqcount'`

- [ ] **Step 3: 写最小实现**

创建 `src/zfc/fastqcount.py`：

```python
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
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (2 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): add reverse_complement and gzip-aware open"
```

---

## Task 2: `_load_library`（跳过 header、两种模式、first-wins）

**Files:**
- Modify: `src/zfc/fastqcount.py`
- Test: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

把测试文件顶部 import 改为：

```python
from zfc.fastqcount import (
    _load_library,
    _open_text,
    reverse_complement,
)
```

在 `test_open_text_detects_gzip` 之后追加：

```python
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
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ImportError: cannot import name '_load_library'`

- [ ] **Step 3: 实现**

在 `_open_text` 之后追加：

```python
def _load_library(path: str, mode: str):
    """Parse the library file into an entry list and a lookup index.

    The first line is skipped when it looks like a header (first field is
    ``gene``). ``index`` maps a guide (guide mode) or a (guide, barcode)
    (barcode mode) to a gene; the first occurrence wins on duplicates.
    """
    entries = []
    index = {}
    with open(path) as handle:
        lines = [line.rstrip("\r\n") for line in handle]
    if not lines:
        return entries, index
    first_col = lines[0].split("\t", 1)[0]
    data = lines if first_col != "gene" else lines[1:]
    for lineno, line in enumerate(data, 1):
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
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): parse library file into entries and index"
```

---

## Task 3: `_iter_paired_records` 与 `_iter_batches`

**Files:**
- Modify: `src/zfc/fastqcount.py`
- Test: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

import 改为：

```python
from zfc.fastqcount import (
    _iter_paired_records,
    _load_library,
    _open_text,
    reverse_complement,
)
```

追加测试：

```python
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
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ImportError: cannot import name '_iter_paired_records'`

- [ ] **Step 3: 实现**

在 `_load_library` 之后追加：

```python
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
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (7 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): stream paired FASTQ records and batches"
```

---

## Task 4: `_match_guide` 与 `_match_barcode_pair`

**Files:**
- Modify: `src/zfc/fastqcount.py`
- Test: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

import 改为：

```python
from zfc.fastqcount import (
    _iter_paired_records,
    _load_library,
    _match_barcode_pair,
    _match_guide,
    _open_text,
    reverse_complement,
)
```

追加测试：

```python
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
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ImportError: cannot import name '_match_guide'`

- [ ] **Step 3: 实现**

在 `_iter_batches` 之后追加：

```python
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
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (10 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): extract guides and barcodes from read pairs"
```

---

## Task 5: worker 与 `_count`（串行路径 + 并行路径）

**Files:**
- Modify: `src/zfc/fastqcount.py`
- Test: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

import 改为：

```python
from zfc.fastqcount import (
    _count,
    _iter_paired_records,
    _load_library,
    _match_barcode_pair,
    _match_guide,
    _open_text,
    reverse_complement,
)
```

追加测试：

```python
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
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ImportError: cannot import name '_count'`

- [ ] **Step 3: 实现**

在 `_match_barcode_pair` 之后追加 worker 全局与 `_count`：

```python
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
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (12 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): count batches serially or via process pool"
```

---

## Task 6: `_write_output`

**Files:**
- Modify: `src/zfc/fastqcount.py`
- Test: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

import 改为：

```python
from zfc.fastqcount import (
    _count,
    _iter_paired_records,
    _load_library,
    _match_barcode_pair,
    _match_guide,
    _open_text,
    _write_output,
    reverse_complement,
)
```

追加测试：

```python
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
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ImportError: cannot import name '_write_output'`

- [ ] **Step 3: 实现**

在 `_count` 之后追加：

```python
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
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (15 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): write headered count tables"
```

---

## Task 7: `_build_patterns`、`main` 与 `__main__`

**Files:**
- Modify: `src/zfc/fastqcount.py`
- Test: `tests/test_fastqcount.py`

- [ ] **Step 1: 写失败测试**

import 改为：

```python
from zfc.fastqcount import (
    _count,
    _iter_paired_records,
    _load_library,
    _match_barcode_pair,
    _match_guide,
    _open_text,
    _write_output,
    main,
    reverse_complement,
)
```

追加测试：

```python
def test_cli_help_exits_zero():
    with pytest.raises(SystemExit) as exc_info:
        main(["--help"])
    assert exc_info.value.code == 0


def test_cli_missing_read_errors(tmp_path):
    with pytest.raises(SystemExit) as exc_info:
        main(["-1", str(tmp_path / "nope.fastq"), "-2", str(tmp_path / "nope2.fastq")])
    assert exc_info.value.code == 2


def test_cli_invalid_pattern_errors(tmp_path):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, ["AAA"])
    _write_fastq(r2, ["TTT"])
    with pytest.raises(SystemExit) as exc_info:
        main(["-1", str(r1), "-2", str(r2), "-p", "ACCGTTT"])
    assert exc_info.value.code == 2


def test_cli_end_to_end_guide_raw(tmp_path, capsys):
    r1 = tmp_path / "r1.fastq"
    r2 = tmp_path / "r2.fastq"
    _write_fastq(r1, [_guide_seq(GUIDE_A), _guide_seq(GUIDE_A)])
    _write_fastq(r2, ["T" * 20, "T" * 20])
    main(["-1", str(r1), "-2", str(r2), "--workers", "1", "--quiet"])
    captured = capsys.readouterr()
    assert captured.out == "guide\tcount\n{0}\t2\n".format(GUIDE_A)
```

- [ ] **Step 2: 运行确认失败**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: FAIL — `ImportError: cannot import name 'main'`

- [ ] **Step 3: 实现**

在 `_write_output` 之后追加：

```python
def _require_groups(pattern, name, minimum):
    compiled = re.compile(pattern)
    if compiled.groups < minimum:
        raise ValueError(
            "{} must define at least {} capture group(s): {!r}".format(
                name, minimum, pattern
            )
        )
    return compiled


def _build_patterns(args):
    if args.mode == "guide":
        _require_groups(args.guide_pattern, "guide pattern", 1)
        return {"guide": args.guide_pattern}
    _require_groups(args.anchor_pattern, "anchor pattern", 1)
    _require_groups(args.barcode_pattern, "barcode pattern", 1)
    return {"anchor": args.anchor_pattern, "barcode": args.barcode_pattern}


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="zfc-count",
        description="Count sgRNAs (and optional iBAR barcodes) from paired-end FASTQ files.",
    )
    parser.add_argument(
        "--mode",
        choices=["guide", "barcode"],
        default="guide",
        help="Counting mode: guide (sgRNA only) or barcode (sgRNA + iBAR).",
    )
    parser.add_argument("-1", "--read1", required=True, help="Forward/read1 FASTQ file.")
    parser.add_argument("-2", "--read2", required=True, help="Reverse/read2 FASTQ file.")
    parser.add_argument(
        "-l", "--library",
        help="Library file (gene<tab>guide, or gene<tab>guide<tab>barcode). "
             "First line may be a header. Omit for raw counts without gene mapping.",
    )
    parser.add_argument("-p", "--guide-pattern", default=DEFAULT_GUIDE, help="sgRNA regex pattern.")
    parser.add_argument("-a", "--anchor-pattern", default=DEFAULT_ANCHOR, help="Full sgRNA+barcode regex pattern.")
    parser.add_argument("-b", "--barcode-pattern", default=DEFAULT_BARCODE, help="Barcode regex pattern.")
    parser.add_argument("-o", "--output", default="-", help="Output file, or '-' for stdout.")
    parser.add_argument("--workers", type=int, default=0, help="Worker processes; 0 = auto, 1 = serial.")
    parser.add_argument("--chunk-size", type=int, default=50000, help="Read pairs per worker batch.")
    parser.add_argument("--quiet", action="store_true", help="Suppress the stderr summary line.")
    args = parser.parse_args(argv)

    for path in (args.read1, args.read2):
        if not os.path.isfile(path):
            parser.error("file not found: {}".format(path))
    if args.library is not None and not os.path.isfile(args.library):
        parser.error("file not found: {}".format(args.library))

    try:
        pattern_strings = _build_patterns(args)
    except (re.error, ValueError) as exc:
        parser.error(str(exc))

    entries = None
    library = None
    if args.library is not None:
        try:
            entries, library = _load_library(args.library, args.mode)
        except ValueError as exc:
            parser.error(str(exc))

    workers = args.workers if args.workers > 0 else min(os.cpu_count() or 1, 8)
    counts, pairs = _count(
        args.mode, pattern_strings, library,
        args.read1, args.read2, workers, args.chunk_size,
    )
    if not args.quiet:
        print("zfc-count: processed {} read pairs".format(pairs), file=sys.stderr)
    _write_output(counts, entries, args.mode, args.output)


if __name__ == "__main__":
    main()
```

- [ ] **Step 4: 运行确认通过**

Run: `uv run pytest tests/test_fastqcount.py -v`
Expected: PASS (19 passed)

- [ ] **Step 5: 提交**

```bash
git add src/zfc/fastqcount.py tests/test_fastqcount.py
git commit -m "feat(fastqcount): add zfc-count CLI"
```

---

## Task 8: 注册 console script + 文档更新 + 全量校验

**Files:**
- Modify: `pyproject.toml`
- Modify: `docs/data-format.md`
- Modify: `docs/architecture.md`
- Modify: `docs/development.md`
- Modify: `README.md`

- [ ] **Step 1: 注册 console script**

`pyproject.toml` 里把：

```toml
[project.scripts]
zfc = "zfc.cli:main"
```

改为：

```toml
[project.scripts]
zfc = "zfc.cli:main"
zfc-count = "zfc.fastqcount:main"
```

- [ ] **Step 2: 更新 `docs/data-format.md`**

把 `## bin/ counting scripts` 一整节替换为：

```markdown
## FASTQ counting (zfc-count)

The legacy bash tools `bin/library_count_sgrna` and
`bin/library_count_sgrna_with_barcode` are deprecated in favor of the Python
`zfc-count` CLI (module `src/zfc/fastqcount.py`):

```
zfc-count --mode {guide,barcode} -1 R1.fastq[.gz] -2 R2.fastq[.gz] [-l library.txt]
```

- Without `-l` → raw counts: `guide<TAB>count` (guide mode) or
  `guide<TAB>barcode<TAB>count` (barcode mode); not directly consumable by
  `zfc` (no `gene` column).
- With `-l` → `gene<TAB>guide<TAB>barcode<TAB>count`, matching the `zfc`
  input shape (counts become ctrl/exp after pivoting samples). Guide mode
  emits `barcode == guide`.
- Library (guide mode): `gene <tab> guide`; (barcode mode):
  `gene <tab> guide <tab> barcode`. A header row beginning with `gene` is
  skipped.
- gzip input is auto-detected by magic bytes (no extension required).
- Patterns default to the iBAR design; overridable via
  `-p/-a/-b`. Fuller details: `docs/specs/2026-08-24-fastqcount-design.md`.
```

- [ ] **Step 3: 更新 `docs/architecture.md`**

把顶部布局里的这一行：

```
- `bin/` — three executable entry points (one Python CLI, two bash scripts).
```

改为：

```
- `bin/` — repo-only tools: two legacy bash counting scripts (deprecated;
  superseded by the `zfc-count` console script in `src/zfc/fastqcount.py`).
```

把 `## bin/ entry points` 一节开头补一句（放在「The Python CLI is the ...」之前，或直接在其后追加一行）：

```
The console script `zfc-count` (module `src/zfc/fastqcount.py`) replaces the
two bash counting scripts below; both remain only for backward reference.
```

把 `## Where to add code` 最后一条：

```
- New top-level command → `bin/` (later registered as a `[project.scripts]` entry).
```

改为：

```
- New top-level command → add it under `src/zfc/` and register it in
  `[project.scripts]` (see `src/zfc/fastqcount.py` + the `zfc-count` entry).
```

- [ ] **Step 4: 更新 `docs/development.md`**

在 `## Run the CLI` 一节，把：

```
- `library_count_sgrna` and `library_count_sgrna_with_barcode` are bash scripts
  run from the repo as `./bin/...`; see `docs/data-format.md` for their I/O.
```

改为：

```
- `uv run zfc-count --help` for FASTQ counting (replaces the deprecated
  `bin/library_count_sgrna*` bash scripts); see `docs/data-format.md` for I/O.
```

- [ ] **Step 5: 更新 `README.md`（可选，但推荐）**

在「Quick Start」之后新增一个短节：

```markdown
### Count reads from FASTQ (zfc-count)

`zfc-count` replaces the legacy `bin/library_count_sgrna[,_with_barcode]` bash
scripts. It streams paired-end FASTQ (plain or gzip) and emits the count table
that `zfc` consumes:

```bash
uv run zfc-count --mode barcode -1 R1.fastq.gz -2 R2.fastq.gz -l library.txt -o counts.tsv
```

Omit `-l` to emit raw `guide`/`barcode` counts without gene mapping. See
`docs/data-format.md`.
```

- [ ] **Step 6: 安装并全量校验**

Run:

```bash
uv sync
uvx ruff check
uv run pytest
uv run zfc-count --help
```

Expected: `uv sync` 成功；`ruff` 报 0 个错误；`pytest` 全绿；`zfc-count --help` 打印 usage 并退出 0。

- [ ] **Step 7: 提交**

```bash
git add pyproject.toml uv.lock docs/data-format.md docs/architecture.md docs/development.md README.md
git commit -m "docs(fastqcount): register zfc-count and deprecate bash counters"
```

> 注：`uv sync` 若改动 `uv.lock` 则一并提交（新脚本无新依赖，一般不会变）。

---

## Task 9: benchmark 脚本 + 真实数据吞吐

**Files:**
- Create: `scripts/bench_fastqcount.py`

- [ ] **Step 1: 写 benchmark 脚本**

```python
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
```

- [ ] **Step 2: 快速 smoke（小样本）**

先生成 100 万 read 的子集（两台 risp 各 4,000,000 行）：

```bash
mkdir -p tmp/sra_tools/subset
gzip -dc tmp/sra_tools/fastq/SRR7975595_pass_1.fastq.gz | head -n 4000000 > tmp/sra_tools/subset/r1.fastq
gzip -dc tmp/sra_tools/fastq/SRR7975595_pass_2.fastq.gz | head -n 4000000 > tmp/sra_tools/subset/r2.fastq
```

跑 smoke：

```bash
uv run python scripts/bench_fastqcount.py -1 tmp/sra_tools/subset/r1.fastq -2 tmp/sra_tools/subset/r2.fastq --mode guide --workers 1
uv run python scripts/bench_fastqcount.py -1 tmp/sra_tools/subset/r1.fastq -2 tmp/sra_tools/subset/r2.fastq --mode guide
```

Expected: 打印 `mode=... pairs=... rate=...`；`--workers 1` 与默认（多核）各自成型，默认应更快。

- [ ] **Step 3: 全量基准（可选，耗时长）**

```bash
uv run python scripts/bench_fastqcount.py \
  -1 tmp/sra_tools/fastq/SRR7975595_pass_1.fastq.gz \
  -2 tmp/sra_tools/fastq/SRR7975595_pass_2.fastq.gz \
  --mode guide
```

Expected: 输出 `pairs≈50_049_729`（以 pass_1 实际为准）与 `rate`，记录作为基线（后续对比旧 bash 管线）。

- [ ] **Step 4: 最终校验 + 提交**

```bash
uvx ruff check
uv run pytest
git add scripts/bench_fastqcount.py
git commit -m "feat(fastqcount): add throughput benchmark script"
```

Expected: ruff 0 错误；pytest 全绿。

---

## Self-Review 记录

- **Spec 覆盖**：guide/barcode 两种模式（Task 4/5/7）、raw/zfc-ready（Task 5/6/7）、gzip（Task 1/3）、并行+serial（Task 5）、统一 CLI（Task 7）、console script 注册 + 文档（Task 8）、测试与基准（Task 1–7 测试、Task 9 基准）、不引依赖（Task 8 `uv sync` 无新依赖）。
- **一致性**：函数名在「文件结构」表、各 Task 实现与测试 import 中一致；`Counter` 键约定（guide 串 / `(gene, guide)` / `(gene, guide, barcode)` / `(guide, barcode)`）在 `_process_batch` 与 `_write_output` 两侧对齐。
- **无占位符**：所有代码块为完整内容，命令含预期输出。