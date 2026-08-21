# Algorithm & Statistics

When to read: before modifying fold-change/statistics code, the null
distribution, p-value/FDR handling, or rank aggregation.

Authoritative references: `README.md` ("Algorithm"), `doc/algorithm.tex`, and
the cited papers (Xu et al. 2021; Zhu et al. 2019; Kolde et al. 2012).

## Pipeline (`src/zfc/zfoldchange.py::zfoldchange`)

1. **Normalization** — `df_normalization` (`src/zfc/statsfunc.py`); small-count
   bump: try quantiles 0.05, 0.10, 0.15 (first non-empty wins after dropping
   zeros), else fall back to the 15% quantile of positive ctrl. The mean small
   count is added to every ctrl/exp count.
2. **Fold change** — `fc = exp / ctrl`, `lfc = log2(fc)`.
3. **LFC std model** — keep rows inside lfc quantiles [0.05, 0.95] and ctrl
   quantiles [0.025, 0.95]; bin ctrl into 200 bins; fit a linear regression of
   per-bin lfc std vs bin ctrl mean; predict `lfc_std` per row; clamp
   non-positive `lfc_std` to the median, and rows above the max training ctrl
   to the median `lfc_std`.
4. **Raw ZLFC** — `zlfc = lfc / lfc_std`; up/down ranks are computed here.
5. **sgRNA level** — mean `zlfc` per (gene, guide); empirical p via a per-count
   ECDF null built from `iteration` shuffles; FDR = BH.
6. **Gene level** — mean `zlfc` per gene; same empirical-p + BH approach.
7. **Rank aggregation** — robust rank aggregation (RRA) and mean-rank
   aggregation per sgRNA and per gene, down and up; RRA scores are BH-adjusted.

Returns `(bar_df, sg_df, g_df, train_data, reg)`.

## Normalization methods (`src/zfc/statsfunc.py::df_normalization`)

`total` (default), `median`, `upper_quartile` (0.75), `median_ratio`, `none`.
Implemented by `df_total_count_normfactor`, `df_median_normfactor`,
`df_quantile_normfactor`, `df_median_ratio_normfactor`.

### `upper_quartile` branch (fixed in 0.2.0)

`df_normalization` previously had a dead `elif method == 'quantile':` branch, so
`--normalization upper_quartile` silently fell through to `total`. Fixed: the
branch now checks `upper_quartile` and routes to `df_quantile_normfactor`
(covered by `tests/test_normalization.py::test_upper_quartile_uses_quantile_factor`).

## Statistical helpers (`src/zfc/statsfunc.py`)

- `ecdf(x)` — empirical CDF returned as a callable (searchsorted over sorted x);
  retained as a helper, though `zfoldchange` now uses the equivalent vectorized
  `np.searchsorted(... , side='left')` directly.
- `p_adjust(p, method)` — `'BH'` (`bh_correction`) or `'Bonferroni'`.
- `df_robust_rank_aggregation(rank_matrix)` — RRA via beta-score rho (used); vectorized.
- `df_mean_rank_aggregation(rank_matrix)` — mean-rank normal approximation (used).

## Defaults & knobs

- `--null-iteration` default `100` (higher = slower, better-calibrated p).
- `--top-n-sgrna` / `--top-n-gene` default `None` → the median per-group barcode
  count is used as the rank-matrix width.

## Resolved issues (0.2.0)

All previously-documented issues are fixed and covered by tests:

1. **Non-reproducible p-values** — fixed by threading a `seed` through
   `zfoldchange` (default `0`; `-1` = nondeterministic) via
   `np.random.default_rng`; `tests/test_zfoldchange.py::test_seed_determinism`.
2. **`median_ratio` failure called `exit()`** — now `warnings.warn(...)` and
   falls back to `total`;
   `tests/test_normalization.py::test_median_ratio_falls_back_to_total_on_many_zeros`.
3. **`df_rank_order_aggregation` used `np.math.factorial`** (removed in NumPy 2)
   — function deleted.
4. **Division by zero unguarded** — input validation (numeric, non-negative,
   at-least-one-positive) plus a post-normalization positive-`ctrl` guard;
   `tests/test_zfoldchange.py::test_zero_counts_do_not_produce_inf`,
   `test_negative_counts_raise`, `test_all_zero_column_raises`.
5. **`upper_quartile` fell through to `total`** — see the `upper_quartile`
   section above.

Change any numeric behavior only with a matching test (`docs/development.md`).