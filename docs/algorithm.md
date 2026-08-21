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

### Known bug: `upper_quartile` branch never runs

- The CLI (`bin/zfc`) and the `df_normalization` assert accept `upper_quartile`,
  but the code branch that applies the quantile factor checks `quantile`
  (`elif method == 'quantile':`). Net effect: `--normalization upper_quartile`
  silently falls through to `total` normalization. `README.md` also still shows
  the old name `quantile`.
  - Why (source): a rename between README/CLI and `statsfunc` was not completed.
  - When: relevant if anyone relies on upper-quartile normalization — fix it
    and add a test.
  - Expires: when the name is unified (prefer `upper_quartile`) and covered by a test.

## Statistical helpers (`src/zfc/statsfunc.py`)

- `ecdf(x)` — empirical CDF returned as a callable (searchsorted over sorted x).
- `p_adjust(p, method)` — `'BH'` (`bh_correction`) or `'Bonferroni'`.
- `df_robust_rank_aggregation(rank_matrix)` — RRA via beta-score rho (used).
- `df_mean_rank_aggregation(rank_matrix)` — mean-rank normal approximation (used).
- `df_rank_order_aggregation(rank_matrix)` — dead code; unused by the pipeline.

## Defaults & knobs

- `--null-iteration` default `100` (higher = slower, better-calibrated p).
- `--top-n-sgrna` / `--top-n-gene` default `None` → the median per-group barcode
  count is used as the rank-matrix width.

## Known issues (fix deliberately; do not paper over)

1. **Non-reproducible p-values.** Null distributions are built with
   `df.sample(frac=1)` and no `random_state`, so p/FDR/RRA can vary between
   runs of the same input. If reproducibility is requested, thread a fixed seed
   through `zfoldchange`.
   - Why: observed in `zfoldchange.py` (two `.sample(frac=1)` calls).
   - When: reproducibility/regression work.
   - Expires: when a seed is threaded through and locked by a test.

2. **`median_ratio` failure calls `exit()`.** Both failure branches terminate
   the process and print a misleading "using total count normalization" message
   even though they do _not_ fall back. They should raise or actually fall back.
   - Why: `src/zfc/statsfunc.py::df_normalization` lines that call `exit(...)`.
   - When: any normalization-quality or error-handling work.
   - Expires: when replaced with raise/fallback and tested.

3. **`df_rank_order_aggregation` uses `np.math.factorial`**, removed in
   NumPy >= 2.0. It is dead code — delete it, or port it only if ever needed.
   - Why: `np.math` is gone in NumPy 2; the function is not imported/used.
   - When: NumPy 2 compatibility or dead-code cleanup.
   - Expires: when the function is deleted or ported.

4. **Division by zero is unguarded.** `fc = exp / ctrl` turns ctrl=0 (raw 0, or
   NaN replaced by 0 up front) into inf/nan in `lfc`; the small-count bump
   mitigates positive counts but there is no explicit guard.
   - Why: `data.replace(np.nan, 0)` then direct division in `zfoldchange.py`.
   - When: robustness/edge-case work.
   - Expires: when an explicit zero-handling policy is added and tested.

Change any numeric behavior only with a matching test (`docs/development.md`).