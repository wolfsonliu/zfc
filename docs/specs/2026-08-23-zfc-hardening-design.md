# Design: ZFC Hardening — Robustness, Code Quality, Engineering (Output-Preserving)

## Problem Statement

`zfc` is correct and, after the 0.2.0 vectorization pass, fast enough. What
remains is a cluster of smaller issues that reduce robustness and
maintainability without touching the statistics:

- Running `zfc` with no `-i` reaches `os.path.isfile(None)` and crashes with a
  `TypeError` traceback instead of a clean usage error.
- `zfoldchange()` validates its required columns with a bare `assert`, which is
  stripped under `python -O` and inconsistent with `read_raw_count`'s
  `ValueError`.
- `zfoldchange()` is a 313-line, seven-step monolith; two of its steps (sgRNA
  and gene empirical p-values) are near-duplicates.
- `plot.py` repeats the same 2×4 scatter/histogram layout four times and has a
  misleadingly named parameter.
- There is no CI, and only 10 tests — nothing guards against a refactor
  silently changing a floating-point result.
- `README.md` documents a stale help with removed flags (`--punish-rate`,
  `--n-sd`).

The hard constraint for this round: **numeric output must not change**.
Column names, file names, sort order, and every numeric value must stay
bit-for-bit identical to `0.2.1`.

## Proposed Solution

Do three pillars of work in increasing risk order, each guarded by a new
"golden parity" test that locks the current numeric output as a regression
baseline. This is an output-preserving refactor plus engineering hardening —
no statistical-convention changes, no performance work, no new features.

## Decisions (locked via Q&A)

| Question | Decision |
| --- | --- |
| Optimization focus | Correctness/robustness + code quality/architecture + engineering/test/CI |
| Feature expansion | **Out of scope this round** |
| Performance/parallelism | **Out of scope this round** |
| Numeric output | **Strictly unchanged** (bit-for-bit); statistical-convention fixes deferred |
| Execution strategy | Plan A — incremental, output-preserving, parity-guarded |
| Version | No release version change this round (output contract unchanged) |

## Detailed Design

### Pillar 1 — Robustness (no output change)

1. **`cli.py` — make `--input` required.** Today `-i/--input` has no `required`
   and no default, so `args['input']` is `None` and `os.path.isfile(None)`
   raises `TypeError`. Set `required=True` on the argument so argparse prints a
   standard "the following arguments are required: -i/--input" error and exits
   with code 2. Keep the existing `os.path.isfile` check for a friendly
   "input file not found" message on a missing-but-provided path.

2. **`zfoldchange.py` — replace the column `assert` with `ValueError`.** The
   loop at the top of `zfoldchange` (`assert x in data.columns, ...`) becomes a
   `ValueError` listing the missing column(s), mirroring `read_raw_count`. The
   identifier/count validation that follows is already `ValueError`-based and
   stays as-is.

3. **Add a golden parity test (the safety net).** `tests/test_parity.py` runs
   `zfoldchange` on a fixed fixture with a fixed `seed`, serializes each of the
   three result tables with a fixed float format
   (`to_csv(float_format='%.17g', sep='\t')`), and asserts the resulting bytes
   equal reference bytes captured once from the current `0.2.1` implementation
   (stored as a small text fixture in the repo). This test must be added
   **before** any refactor so it can prove the refactor is output-neutral. It
   becomes the executable form of "strictly no numeric change."

### Pillar 2 — Code quality / architecture (behavior-identical refactor)

4. **Split `zfoldchange()` into private step functions** (same module), keeping
   the public signature and the returned `(bar_df, sg_df, g_df, train_data,
   reg)` tuple unchanged:

   - `_validate_and_normalize(data, normalization)` → `(bar_df, smallcount)`
   - `_fold_change(bar_df)`
   - `_lfc_std_model(bar_df)` → `(bar_df, train_data, reg)`
   - `_barcode_zlfc_rank(bar_df)` (adds `zlfc`, `rank_down`, `rank_up`, `p`)
   - `_level_pvalue(z, codes, counts, iteration, rng)` — **shared** by the
     sgRNA and gene steps, collapsing the two near-identical permutation/ECDF
     blocks (removes ~40 duplicated lines). Returns the empirical `p` array.
   - `_rank_aggregation(bar_df, level, top_n)` — returns the `RRA_*`/`Mean_Rank_*`
     columns for a level, collapsing the sgRNA and gene rank-matrix blocks.

   Each helper moves code as-is (no reordering of float operations), so the
   parity test stays green. `mean()` and the `bincount` null loop are preserved
   exactly; the RNG draw order (`rng.permutation` once per iteration) is
   preserved too.

5. **Type annotations.** Add `from __future__ import annotations` and annotate
   the public API (`zfoldchange`, `read_raw_count`, `write_bar_result`,
   `write_sg_result`, `write_g_result`) and the statsfunc helpers. Uses
   `X | None` (Python >= 3.11, per `pyproject.toml`).

6. **Docstrings.** Give math helpers a one-line reference citation where
   missing, matching the existing pattern in `src/zfc/statsfunc.py`.

7. **`plot.py` cleanup.**
   - Extract the shared 2×4 scatter/histogram layout used by
     `normcount_scatter`, `lfc_normcount_scatter`, `zlfc_normcount_scatter`
     into one parameterized builder (da/dx/dy + labels).
   - Rename `counts_boxplot(inputdata, sgresult)`'s `sgresult` parameter to
     `barresult` (it is called with `barresult`, the normalized counts frame).
   - Vectorize `gini` (`np.dot` instead of the Python `sum` loop). Note: the
     value feeds only a plot title rounded to 4 decimals, so it never enters a
     result table; a last-ulp float-order difference is invisible.

### Pillar 3 — Engineering / test / CI

8. **Add `.github/workflows/ci.yml`** — matrix over Python 3.11/3.12/3.13,
   steps: install uv, `uv sync`, `uvx ruff check`, `uv run pytest`.

9. **Expand tests** (`tests/test_cli.py`, `tests/test_zfoldchange.py`,
   `tests/test_parity.py`):
   - `main([])` (no `-i`) raises `SystemExit` code 2 (no `TypeError`).
   - `zfoldchange` missing a required column raises `ValueError` (not
     `AssertionError`).
   - Golden parity: refactored output equals the 0.2.1 baseline.

10. **Fix `README.md` drift.** Remove the stale `--punish-rate`/`--n-sd` from
    the help block and re-sync it to the current argparse output; update the
    Dependency section to uv + `requires-python >= 3.11`.

## Out of Scope (this round)

- **Performance**: `--workers`/multiprocessing, numba/JIT, further
  vectorization of the null loop or rank-matrix build.
- **Statistical-convention changes** (numeric output would change; deferred to
  a separate, author-confirmed design):
  - the two-sided p reflection `min(p, 1-p)` (half of the conventional
    two-sided p) at barcode/sgRNA/gene levels;
  - the permutation ECDF `side='left'` convention (`P(null < x)` vs
    `P(null ≤ x)`);
  - `lfc_std` clamped to the median (vs. linear extrapolation) for ctrl values
    above the training range.
- **Feature expansion**: new normalization methods, new output columns/files,
  new plots, and DESeq2-style multi-replicate variance modeling.

## Success Criteria

- `uvx ruff check` clean; `uv run pytest` green (existing + new tests).
- Golden parity test passes: barcode/sgRNA/gene tables are bit-for-bit
  identical to the 0.2.1 baseline after the Pillar 2 refactor.
- `zfc` with no `-i` prints a clean argparse error and exits code 2 (no
  traceback).
- CI runs lint + tests on push/PR for Python 3.11–3.13.
- No change to any output filename, column, sort order, or numeric value.

## Implementation Order

1. Pillar 1 (robustness) + the golden parity test, so the baseline is locked
   before anything else moves.
2. Pillar 2 (refactor), re-running the parity test after each helper split.
3. Pillar 3 (CI, tests, README).

## Open Questions

- None blocking for the output-preserving scope. The three deferred
  statistical-convention items (p reflection, ECDF side, `lfc_std`
  extrapolation) are recorded above for a future, author-reviewed design; they
  are intentionally not touched here.