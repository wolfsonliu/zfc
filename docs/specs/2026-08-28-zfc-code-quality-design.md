# Design: ZFC Code Quality & Robustness Hardening (Output-Preserving)

## Problem Statement

`zfc` is correct and, after the 0.2.0 vectorization and 0.2.1 hardening
rounds, fast enough (`zfoldchange(iteration=100)` ≈ 3.77s on
`tmp/REP_A549_High.txt`, ~233k barcodes). What remains is a small cluster of
issues that reduce maintainability and robustness without touching the
statistics:

- `statsfunc.py::ecdf` is dead code (the pipeline uses the equivalent
  vectorized `np.searchsorted(..., side='left')` directly).
- `p_adjust` validates its `method` argument with a bare `assert`, which is
  stripped under `python -O` and inconsistent with the `ValueError` policy
  used elsewhere in the package.
- Type annotations are incomplete: `zfoldchange`, `_lfc_std_model`, and
  `_empirical_p` lack return types.
- Two plot defects: `zlfc_rra_scatter` / `zlfc_rrafdr_scatter` compute
  `np.log10(rra) * -1` and emit `inf` + a `RuntimeWarning` when an RRA score
  is exactly 0; `lfcstd_ctrlmean_figure` uses `np.arange(min, max)`, which
  degenerates to a single-point line when the control-mean range spans less
  than one unit.
- `_lfc_std_model` runs `pd.cut(...).astype(str)` over every row (a
  categorical→string conversion on ~194k elements ≈ 0.55s) solely to get
  string bin labels for `groupby`.

Hard constraint for this round: **numeric output must not change.** Column
names, file names, sort order, and every numeric value of the barcode/sgRNA/
gene result tables must stay bit-for-bit identical to `0.2.2`. The existing
golden parity test (`tests/test_parity.py`) is the safety net.

## Proposed Solution

Six small, output-preserving changes to `statsfunc.py`, `zfoldchange.py`, and
`plot.py`, plus two focused tests. No statistical-convention changes, no new
features, no new runtime dependencies, no change to the public output contract.

## Decisions (locked via Q&A)

| Question | Decision |
| --- | --- |
| Optimization focus | Code quality / robustness (safe, zero-risk) |
| Item scope | All six items in the candidate list (incl. the safe `pd.cut` micro-optimization) |
| Numeric output | Strictly unchanged (bit-for-bit), parity-guarded |
| Performance work | Only the incidental `pd.cut` win; no `--workers`/numba/rank-agg rewrite |
| Statistical-convention fixes | Out of scope (still deferred) |

## Detailed Design

### 1. `statsfunc.py` — delete dead `ecdf`

Remove the `ecdf` function (lines 16–36) and its now-unused `import functools`
(line 11). No callers exist in `src/` or `tests/` (verified by grep); the
pipeline computes empirical CDF values inline at `zfoldchange.py:183`.

### 2. `statsfunc.py` — `p_adjust` `assert` → `ValueError`

Replace

```python
assert method in ['BH', 'Bonferroni'], 'Adjust method should be in: BH, Bonferroni'
```

with

```python
if method not in ('BH', 'Bonferroni'):
    raise ValueError('Adjust method should be in: BH, Bonferroni')
```

Callers (`zfoldchange.py`) pass only `'BH'`; behavior is unchanged. Survives
`python -O` and matches the package's `ValueError`-based validation idiom.

### 3. `zfoldchange.py` — complete type annotations

Add return annotations (no new imports; `sklearn` is already imported):

- `zfoldchange(...) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, linear_model.LinearRegression]`
- `_lfc_std_model(bar_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, linear_model.LinearRegression]`
- `_empirical_p(groups, values: pd.Series, iteration: int, rng) -> pd.DataFrame`
- `_rank_aggregation(...) -> dict[str, pd.Series]`

Leave the `groups` argument of `_empirical_p` unannotated: there is no clean,
stable public `SeriesGroupBy` type across the supported pandas versions.

### 4. `plot.py` — RRA scatter uses `_log10p`

In `zlfc_rra_scatter` and `zlfc_rrafdr_scatter`, replace
`np.log10(rra) * -1` with the existing `_log10p(rra)` helper. For `rra == 0`
this yields the large-but-finite `-log10(tiny)` top value instead of `inf` and
suppresses the `RuntimeWarning`; for `rra > 0` the value is unchanged. Only
the plot y-values change (never a result-table value).

### 5. `plot.py` — `lfcstd_ctrlmean_figure` uses `np.linspace`

Replace `x = np.arange(train_data['ctrlmean'].min(), train_data['ctrlmean'].max())`
with

```python
x = np.linspace(train_data['ctrlmean'].min(), train_data['ctrlmean'].max(), 100)
```

so the regression line is always rendered as a smooth 100-point line even when
the control-mean range is narrower than one unit.

### 6. `zfoldchange.py` — group on `pd.cut` categorical instead of `.astype(str)`

Today:

```python
bins = pd.cut(model_data['ctrl'], 200).astype(str)
train_data = pd.DataFrame({
    'lfcstd': model_data.groupby(bins)['lfc'].std(),
    'ctrlmean': model_data.groupby(bins)['ctrl'].mean(),
})
```

The `.astype(str)` runs an expensive elementwise categorical→string conversion
on ~194k rows only to produce group keys. Instead, group on the Categorical
and stringify the resulting *index* (≈200 labels):

```python
bins = pd.cut(model_data['ctrl'], 200)
train_data = pd.DataFrame({
    'lfcstd': model_data.groupby(bins, observed=True)['lfc'].std(),
    'ctrlmean': model_data.groupby(bins, observed=True)['ctrl'].mean(),
})
train_data.index = train_data.index.astype(str)
```

`model_data['ctrl']` is already bounded by the [0.025, 0.95] quantiles before
cutting, so the cut has no out-of-range values. Verified on
`tmp/REP_A549_High.txt`: both paths yield the same 200 bins, the same index,
and 0 NaN cuts. A fixed-fixture equivalence test asserts the new `train_data`
(index + values) equals the pre-change result, and the golden parity test
guards the downstream `lfc_std`/`zlfc` chain, so any pathological drift (e.g. a
cut boundary producing NaN on some input) is caught.

## Tests

- `tests/test_statsfunc.py` (new): `p_adjust(..., 'BH')` and
  `'Bonferroni'` return unchanged; `p_adjust(..., 'bogus')` raises `ValueError`.
- `tests/test_zfoldchange.py` (new): `_lfc_std_model` produces a `train_data`
  equivalent to the pre-change `.astype(str)` result on a fixed fixture
  (index labels and both columns equal).
- Existing `tests/test_parity.py` and `tests/test_plot.py` must stay green
  (parity locks the three result tables bit-for-bit; the plot smoke test
  exercises the changed figure builders end-to-end).

## Success Criteria

- `uvx ruff check` clean; `uv run pytest` green (existing + new tests).
- Golden parity test passes: barcode/sgRNA/gene tables bit-for-bit identical
  to the `0.2.2` baseline.
- `train_data` equivalence test passes (item 6 introduces no drift).
- No change to any output filename, column, sort order, or numeric value.

## Out of Scope

- Statistical-convention corrections (two-sided p reflection `min(p, 1-p)`,
  ECDF `side='left'`, `lfc_std` median clamp vs. linear extrapolation) — still
  deferred to a separate, author-confirmed design.
- Performance: `--workers`/multiprocessing, numba/JIT, and full vectorization
  of `_rank_aggregation`.
- New normalization methods, new result columns/files, or new plot types.

## Implementation Order

1. `statsfunc.py`: delete `ecdf` + `functools`; `p_adjust` `ValueError`.
2. `zfoldchange.py`: type annotations; categorical `pd.cut` grouping.
3. `plot.py`: `_log10p` in the two RRA scatters; `np.linspace` in
   `lfcstd_ctrlmean_figure`.
4. Add the two tests; run `ruff` + `pytest` (parity + plot must stay green).

## Open Questions

- None.