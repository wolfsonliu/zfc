# Design: Optimize ZFC — Speed, Efficiency, Bug Fixes

## Problem Statement

`zfc` is correct but slow and has several latent bugs. On real iBAR screening
data (`tmp/REP_A549_{High,HL,Low}.txt`, each ~233k barcode rows → ~58k sgRNAs →
~19k genes), `zfoldchange(..., iteration=100)` extrapolates to **~11 min per
file** (~33 min for all three). Most of that time is spent in Python-level
row-wise loops, and several documented issues (`docs/algorithm.md` "Known
issues") plus a few undocumented inefficiencies reduce correctness and
robustness.

Goal: make the pipeline materially faster, fix the bugs, and keep the public
output contract (columns, filenames, sorting, header names) unchanged.

## Decisions (locked via Q&A)

| Question | Decision |
| --- | --- |
| Priority | Speed + bug fixes together; numeric outputs **may change** (more correct / reproducible) |
| Output contract | Keep columns, filenames, sorting, and `gene/guide/barcode/ctrl/exp` header; only values change |
| Reproducibility | Add `--seed`; default fixed (`0`); `--seed -1` = non-deterministic |
| Approach | A — incremental vectorization (no numba, no multiprocessing this round) |
| Version | `0.1.7` → `0.2.0` |

## Performance Profile (measured on `REP_A549_High.txt`, current code)

Toolchain at measurement time: numpy 2.5.2, scipy 1.18.0, pandas 3.0.5,
scikit-learn 1.9.0. Full run at `iteration=3` = **19.8s**.

| Step | Measured cost | Weight |
| --- | --- | --- |
| RRA rank aggregation (`df_robust_rank_aggregation`, row-wise `apply`) | sg 58243×4 = 6.08s; gene 19363×12 ≈ 2.07s (×2 matrices each ≈ 16s) | **dominant** |
| Null-distribution loop (sgRNA + gene) | 0.12s + 0.02s per iter → ~14s @100 | 2nd |
| LFC-std linear model + 200-bin cut | ~0.59s | minor |
| Normalization + fold change | ~0.04s | negligible |
| `read_raw_count` string index (`apply('_'.join, axis=1)`) | ~0.5s (233k rows) | minor |

## Proposed Solution

Vectorize the two hot paths (RRA and the null-distribution loop) with
numpy/scipy, fix the documented bugs, add the reproducibility `--seed` knob,
and clean up the plot/IO inefficiencies. No new runtime dependencies; the math
for the *unchanged* parts stays numerically identical (locked by tests).

## Detailed Design

### 1. `statsfunc.py` — RRA vectorization

Replace `df_robust_rank_aggregation(rank_matrix)`'s
`rank_matrix.apply(lambda r: rho_score(r), axis=1)` with an all-at-once
computation, keeping the exact same mathematics:

1. `vals = np.asarray(rank_matrix, dtype=float)`; `sorted_vals = np.sort(vals, axis=1)`
   (NaN sorts to the end).
2. `k = (~np.isnan(sorted_vals)).sum(axis=1)` — non-NaN count per row.
3. Build per-column `a = arange(1, p+1)` and per-row/column
   `b = k[:, None] - arange(p)`; `clip(a, 1, None)` and `clip(b, 1, None)` so
   padded (NaN) cells never hit `b <= 0`.
4. `c = scipy.stats.beta.cdf(sorted_vals, a, b)` (one vectorized call).
5. `masked = np.where(np.isnan(sorted_vals), np.nan, c)`;
   `rho = np.nanmin(k[:, None] * masked, axis=1).clip(0, 1)`.
6. Return a `pd.Series(rho, index=rank_matrix.index)`.

This is the same formula as today (`rho = min over j of (k * beta.cdf(r_j, j,
k-j+1))`, clip [0,1]) and must match the old row-wise result **bit-for-bit**;
lock it with a parity test below. (Current code already skips NaN via the same
`~np.isnan` select and `k = sum(~np.isnan(r))`, so behavior is preserved.)

`df_mean_rank_aggregation` is already vectorized; leave it, but reuse its
not-NaN counting idiom for consistency if convenient.

### 2. `statsfunc.py` — bug fixes

- **`upper_quartile` branch** (`df_normalization`): change the dead
  `elif method == 'quantile':` to `elif method == 'upper_quartile':` so the
  `df_quantile_normfactor` path actually runs. The function's assert already
  rejects `'quantile'`, so no aliasing change is needed; this only makes the
  `upper_quartile` name (already accepted by the assert and the CLI) take its
  intended path.
- **`median_ratio` no-exit**: replace both `exit(...)` calls with a real
  fallback to `total` normalization plus `warnings.warn(...)` that the fallback
  happened (message no longer claims it actually falls back when it exits).
  The "too many zeros" and "median factor is zero" branches both warn + fall
  back to total.
- **Delete `df_rank_order_aggregation`** (dead code using `np.math.factorial`,
  removed in NumPy 2). No callers import it.

### 3. `zfoldchange.py` — seed + null loop + guards

- **Seed threading**: add a `seed=0` keyword to `zfoldchange`. Construct
  `rng = np.random.default_rng(seed)` when `seed >= 0`, else
  `np.random.default_rng()` when `seed < 0`. Replace both `.sample(frac=1)`
  calls with `rng.permutation(values)`.
- **Null loop restructure** (sgRNA step, then mirror for gene step):
  - Precompute once, outside the loop: the group-size values
    (`counts_unique`), the grouper mapping row→group, and the observed group
    means/counts (`sg_df`).
  - Per iteration: `perm = rng.permutation(bar_df['zlfc'].to_numpy())` once;
    compute per-group means over the fixed group partition (pandas groupby on a
    precomputed Series, or numpy `bincount`/`add.at` on group ids — the
    group-size map is reused, not rebuilt).
  - Accumulate per-size null arrays; **build each count's ECDF once after the
    loop**, not inside it (today `sg_ecdf_list[c] = ecdf(...)` re-sorts the
    whole accumulated list every iteration).
  - Preserve ECDF semantics exactly: `ecdf(sample) == searchsorted(sorted_x,
    sample, side='left') / n` — keep `side='left'` so the statistical
    convention is unchanged (only the seed and accumulation change).
  - Evaluate observed p via numpy per unique count, then `np.minimum(p, 1-p)`,
    then `p_adjust(..., 'BH')` (unchanged).
- **p-reflection**: replace `p.map(lambda x: x if x <= 0.5 else 1 - x)` (three
  sites) with `np.minimum(p, 1 - p)` — numerically identical.
- **Division-by-zero guard**: after the small-count bump, assert/ensure `ctrl`
  is strictly positive; add input validation (below) that raises a clear error
  on degenerate columns (all-zero / all-NaN) instead of producing `inf` in
  `lfc`.
- **NaN handling**: keep documented `NaN → 0` behavior for `ctrl`/`exp`
  (switch to `DataFrame.fillna(0)` to avoid pandas 3 `replace` surprises), but
  raise if an identifier column (`gene`/`guide`/`barcode`) contains NaN.
- **Input validation** (library and CLI path alike): reject non-numeric or
  negative `ctrl`/`exp`; require at least one positive count per column so
  normalization is defined.

### 4. `programio.py`

- Keep `read_raw_count`'s signature and the `gene_guide_barcode` string index
  (it is asserted in `tests/test_io.py`), but build the index with vectorized
  string concat (`data['gene'] + '_' + data['guide'] + '_' + data['barcode']`)
  instead of `apply('_'.join, axis=1)`.
- `read_raw_count` raises a clear `ValueError` listing missing required columns
  instead of a bare `KeyError`.

### 5. `cli.py`

- Add `--seed` (`type=int`, default `0`; `-1` = random) and pass it through to
  `zfoldchange`.
- Give a friendly error if `--input` does not exist or is unreadable (before
  `read_raw_count`), and keep the "`--outprefix` directory must already exist"
  behavior unchanged (hard constraint #8).
- Update `--help` text for `--normalization` (unchanged names) and the new
  `--seed`.

### 6. `plot.py` — memory / correctness

- Replace `c=['black'] * len(...)` and `edgecolor=['none'] * len(...)` with the
  scalar `c='black'` / `edgecolor='none'` (large per-point lists are wasteful).
- `plt.close(fig)` after each save in `zfc_plot` (14 figures × 2 formats
  currently accumulate).
- In `zlfc_rra_scatter` / `zlfc_rrafdr_scatter`, stop mutating the caller's
  DataFrame via `data.loc[:, 'RRA'] = ...`; compute on a local copy.
- Keep `matplotlib.use('Agg')` (hard constraint #3).

### 7. `__init__.py`

- Bump `__version__` to `0.2.0` (single-sourced; hard constraint #4).

## Tests (add alongside changes)

- `tests/test_statsfunc.py` (new):
  - `upper_quartile` returns `df_quantile_normfactor` output and differs from
    `total` on a non-trivial fixture.
  - `median_ratio` with too many zeros warns and falls back to `total`
    (result equals total; no `SystemExit`).
  - RRA vectorized == old row-wise reference (bit-for-bit) on a fixture with
    uneven NaN counts per row.
- `tests/test_zfoldchange.py`:
  - same `seed` → two runs identical; different `seed` → allowed to differ;
    `seed=-1` runs without error.
  - zero-count robustness: `ctrl=0` rows do not yield `inf` `lfc`; degenerate
    (all-zero) input raises a clear error.
- `tests/test_cli.py`:
  - `--seed 1` accepted (CLI smoke); missing input file errors cleanly.
- `tests/test_io.py`: existing index test still passes (vectorized join is
  string-identical).

## Success Criteria

- `uvx ruff check` clean; `uv run pytest` green (existing + new tests).
- Single file at `iteration=100` drops from ~11 min to **~1–2 min (≥5×)**;
  record before/after timings on `REP_A549_High.txt`.
- RRA output matches old row-wise implementation bit-for-bit (parity test).
- Same `--seed` → identical output twice; columns/filenames/sorting unchanged.
- All five documented bugs (upper_quartile, median_ratio exit, dead
  `np.math.factorial`, non-reproducibility, div-by-zero) addressed and tested.

## Out of Scope (this round)

- Multiprocessing / `--workers` / per-file parallelism (follow-up, after the
  serial hot paths are gone).
- numba JIT or a full numpy rewrite of the pipeline (approach B).
- New normalization methods, new output columns/files, new plot types.
- Interactive plotting (keep headless Agg).
- Fixing the ECDF `side='left'` convention or other statistical-notation
  changes beyond the five named bugs.

## Implementation Order

1. `statsfunc.py`: RRA vectorization + parity test; upper_quartile fix;
   median_ratio fallback; delete `df_rank_order_aggregation`.
2. `zfoldchange.py`: `seed` threading, null-loop restructure, p-reflection,
   div-zero guard + input validation.
3. `programio.py` + `cli.py`: vectorized index, `--seed`, clearer errors.
4. `plot.py`: scalar colors, `plt.close`, no caller mutation.
5. `__init__.py` version bump; update `docs/algorithm.md` (resolve the
   "Known issues" now fixed) and `README.md` normalization wording.
6. `uvx ruff check` + `uv run pytest`; benchmark before/after on
   `REP_A549_High.txt`.

## Open Questions

- None blocking. (`--seed` default fixed at `0` per Q&A; `require a >= 0` for
  counts follows from the non-negativity validation and is covered by tests.)