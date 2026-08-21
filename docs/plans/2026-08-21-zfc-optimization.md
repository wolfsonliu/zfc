# ZFC Optimization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `zfc` ~15–20× faster on real iBAR screening data (≥233k rows; current baseline **82.2s** at `iteration=100` per file), fix five documented bugs, and add deterministic `--seed`, while keeping the output contract (columns, filenames, sorting, header names) unchanged.

**Architecture:** Vectorize the two hot paths in `statsfunc.py` (RRA) and `zfoldchange.py` (permutation null distribution) with numpy/scipy, thread a `seed` through `zfoldchange`, add input validation, fix the two normalization bugs, and clean up plot/IO inefficiencies. No new runtime dependencies; the vectorized RRA is numerically identical to the old row-wise implementation (verified bit-for-bit), and all other math changes are documented bugfixes.

**Tech Stack:** Python 3.11, numpy 2.x, scipy 1.x, pandas 3.x, scikit-learn, matplotlib (Agg). Tests: pytest. Lint: ruff.

**Spec:** `docs/specs/2026-08-21-zfc-optimization-design.md`.

**Note on toolchain:** AGENTS.md prescribes `uv run pytest` / `uvx ruff check`. If `uv`/`uvx` fail with a read-only-cache error in a sandbox, the equivalent fallback is `.venv/bin/python -m pytest` and `.venv/bin/ruff check` (dev deps already in `.venv`).

---

## File Map

| File | Change |
| --- | --- |
| `src/zfc/statsfunc.py` | vectorize RRA; fix `upper_quartile`; `median_ratio` fallback; delete dead code |
| `src/zfc/zfoldchange.py` | seed, vectorized null loop, p-reflection, validation, div-zero guard |
| `src/zfc/programio.py` | vectorized index; clear missing-column error |
| `src/zfc/cli.py` | `--seed`; input-existence error |
| `src/zfc/plot.py` | scalar colors; `plt.close`; stop mutating callers |
| `src/zfc/__init__.py` | version 0.1.7 → 0.2.0 |
| `tests/test_statsfunc.py` | NEW — RRA parity |
| `tests/test_normalization.py` | add upper_quartile / median_ratio tests |
| `tests/test_zfoldchange.py` | add seed/validation/zero tests |
| `tests/test_cli.py` | add `--seed` + missing-input tests |
| `tests/test_plot.py` | NEW — plot smoke test |
| `README.md`, `docs/algorithm.md` | doc updates |

---

### Task 1: Vectorize RRA in `statsfunc.py`

**Files:**
- Create: `tests/test_statsfunc.py`
- Modify: `src/zfc/statsfunc.py:203-220`

- [ ] **Step 1: Write the failing parity test**

Create `tests/test_statsfunc.py`:

```python
import numpy as np
import pandas as pd
from scipy.stats import beta

from zfc.statsfunc import df_robust_rank_aggregation


def _reference_rra(rank_matrix):
    """Independent, straightforward implementation of the documented RRA rho.

    rho = min over j of clip(k * beta.cdf(r_j, j, k-j+1), 0, 1), where r_j are
    the row's sorted non-NaN ranks and k is the number of non-NaN ranks.
    """
    out = []
    for _, row in rank_matrix.iterrows():
        vals = np.sort(row.to_numpy())
        vals = vals[~np.isnan(vals)]
        k = len(vals)
        scores = beta.cdf(vals, np.arange(1, k + 1), np.arange(k, 0, -1))
        out.append((scores * k).clip(0, 1).min())
    return np.array(out)


def test_robust_rank_aggregation_matches_reference():
    rng = np.random.default_rng(0)
    m = pd.DataFrame(rng.uniform(0, 1, (500, 5)))
    m.iloc[::3, 3] = np.nan  # uneven NaN pattern per row
    m.iloc[::5, 4] = np.nan

    result = df_robust_rank_aggregation(m)
    expected = _reference_rra(m)

    assert list(result.index) == list(m.index)
    assert np.array_equal(result.to_numpy(), expected)
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_statsfunc.py::test_robust_rank_aggregation_matches_reference -v`
Expected: PASS (the old implementation is also correct, so this is a *characterization* test that locks behavior). If it does not pass, stop — the reference diverges and must be reconciled before the vectorization.

- [ ] **Step 3: Replace the function with the vectorized form**

In `src/zfc/statsfunc.py`, replace the entire `df_robust_rank_aggregation` body (lines 203–220) with:

```python
def df_robust_rank_aggregation(rank_matrix):
    # Kolde, R., Laur, S., Adler, P., and Vilo, J. (2012). Robust rank
    # aggregation for gene list integration and
    # meta-analysis. Bioinformatics 28, 573–580.
    vals = np.asarray(rank_matrix, dtype=float)
    sorted_vals = np.sort(vals, axis=1)          # NaN sorts to the end
    k = (~np.isnan(sorted_vals)).sum(axis=1)     # non-NaN count per row
    n_cols = sorted_vals.shape[1]

    a = np.clip(np.arange(1, n_cols + 1), 1, None)      # 1..p per column
    b = np.clip(k[:, None] - np.arange(n_cols), 1, None)  # k-j per row/col
    beta_cdf = stats.beta.cdf(sorted_vals, a, b)

    masked = np.where(np.isnan(sorted_vals), np.nan, beta_cdf)
    scores = np.nanmin(np.clip(k[:, None] * masked, 0, 1), axis=1)
    return pd.Series(scores, index=rank_matrix.index)
```

(No change to imports: `pd`, `np`, `stats` are already imported at the top of `statsfunc.py`.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run pytest tests/test_statsfunc.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/test_statsfunc.py src/zfc/statsfunc.py
git commit -m "perf: vectorize robust rank aggregation"
```

---

### Task 2: Fix `upper_quartile` normalization

**Files:**
- Modify: `tests/test_normalization.py`
- Modify: `src/zfc/statsfunc.py:132-133`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_normalization.py` (add `df_quantile_normfactor` to its import line):

```python
from zfc.statsfunc import df_normalization, df_quantile_normfactor
```

```python
def test_upper_quartile_uses_quantile_factor():
    df = pd.DataFrame({
        'ctrl': [10.0, 20.0, 30.0, 25.0],
        'exp': [5.0, 15.0, 25.0, 10.0],
    })
    result = df_normalization(df, 'upper_quartile')
    expected = df.mul(df_quantile_normfactor(df), axis=1)
    assert result.equals(expected)
    assert not result.equals(df_normalization(df, 'total'))
```

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_normalization.py::test_upper_quartile_uses_quantile_factor -v`
Expected: FAIL (current branch checks `quantile`, so `upper_quartile` falls through to `total`; `result.equals(expected)` is False).

- [ ] **Step 3: Fix the branch**

In `src/zfc/statsfunc.py`, change line 132:

```python
    elif method == 'quantile':
```
to
```python
    elif method == 'upper_quartile':
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run pytest tests/test_normalization.py -v`
Expected: PASS (all existing normalization tests too).

- [ ] **Step 5: Commit**

```bash
git add src/zfc/statsfunc.py tests/test_normalization.py
git commit -m "fix: route upper_quartile normalization to quantile factor"
```

---

### Task 3: `median_ratio` fallback instead of `exit()`

**Files:**
- Modify: `tests/test_normalization.py`
- Modify: `src/zfc/statsfunc.py:1-9,124-131`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_normalization.py`:

```python
def test_median_ratio_falls_back_to_total_on_many_zeros():
    df = pd.DataFrame({
        'ctrl': [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        'exp': [0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    })
    with pytest.warns(UserWarning, match='Too many zeros'):
        result = df_normalization(df, 'median_ratio')
    assert result.equals(df_normalization(df, 'total'))
```

(Add `import pytest` at the top of `tests/test_normalization.py` if not present.)

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_normalization.py::test_median_ratio_falls_back_to_total_on_many_zeros -v`
Expected: FAIL (current code calls `exit()`, raising `SystemExit`, not a `UserWarning`).

- [ ] **Step 3: Implement fallback**

Add `import warnings` at the top of `src/zfc/statsfunc.py` (after `import scipy.stats as stats`).

Replace the `median_ratio` branch (lines 124–131):

```python
    elif method == 'median_ratio':
        medianfactor = df_median_ratio_normfactor(df)
        if (medianfactor == 0).any():
            warnings.warn(
                'Median factor is zero; using total count normalization'
            )
        elif ((df == 0).sum(axis=0) / df.shape[0] > 0.45).any():
            warnings.warn(
                'Too many zeros in counts; using total count normalization'
            )
        else:
            normfactor = medianfactor
```

(`normfactor` is already `df_total_count_normfactor(df)` from the top of the function, so warning + falling through leaves it as `total`.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `uv run pytest tests/test_normalization.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/zfc/statsfunc.py tests/test_normalization.py
git commit -m "fix: fall back to total (not exit) on median_ratio failure"
```

---

### Task 4: Delete dead `df_rank_order_aggregation`

**Files:**
- Modify: `src/zfc/statsfunc.py:171-200`

- [ ] **Step 1: Confirm no importers**

Run: `grep -R "df_rank_order_aggregation" src tests`
Expected: only the definition in `src/zfc/statsfunc.py` (no importers).

- [ ] **Step 2: Delete the function**

Delete the entire `df_rank_order_aggregation` function (including its docstring and the `Q(r)` closure) from `src/zfc/statsfunc.py`.

- [ ] **Step 3: Run the suite to confirm still green**

Run: `uv run pytest`
Expected: PASS (nothing imported the removed function).

- [ ] **Step 4: Run lint**

Run: `uvx ruff check src/zfc/statsfunc.py`
Expected: clean.

- [ ] **Step 5: Commit**

```bash
git add src/zfc/statsfunc.py
git commit -m "chore: remove dead df_rank_order_aggregation (np.math removed in numpy 2)"
```

---

### Task 5: Rewrite `zfoldchange` (seed + validation + vectorized null + p-reflection)

**Files:**
- Modify: `tests/test_zfoldchange.py`
- Modify: `src/zfc/zfoldchange.py` (full replacement)

- [ ] **Step 1: Write the failing tests**

Append to `tests/test_zfoldchange.py` (add `import pytest` at the top):

```python
def test_seed_determinism():
    data = _make_input()
    _, sg_a, g_a, _, _ = zfoldchange(data, iteration=5, seed=42)
    _, sg_b, g_b, _, _ = zfoldchange(data, iteration=5, seed=42)
    assert sg_a['p'].equals(sg_b['p'])
    assert g_a['p'].equals(g_b['p'])


def test_negative_seed_uses_entropy():
    data = _make_input()
    zfoldchange(data, iteration=5, seed=-1)  # must not raise


def test_zero_counts_do_not_produce_inf():
    data = _make_input()
    data.loc[0, 'ctrl'] = 0.0
    data.loc[1, 'exp'] = 0.0
    bar, sg, _, _, _ = zfoldchange(data, iteration=5, seed=0)
    assert not np.isinf(bar['lfc']).any()
    assert not np.isinf(bar['zlfc']).any()
    assert not np.isinf(sg['zlfc']).any()


def test_all_zero_column_raises():
    data = _make_input()
    data['ctrl'] = 0.0
    with pytest.raises(ValueError):
        zfoldchange(data, iteration=5, seed=0)


def test_negative_counts_raise():
    data = _make_input()
    data.loc[0, 'ctrl'] = -5.0
    with pytest.raises(ValueError):
        zfoldchange(data, iteration=5, seed=0)
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `uv run pytest tests/test_zfoldchange.py -v`
Expected: `test_seed_determinism` FAIL (TypeError: unexpected keyword `seed`); `test_negative_seed_uses_entropy` FAIL; `test_zero_counts_do_not_produce_inf` FAIL (inf present); `test_all_zero_column_raises` FAIL; `test_negative_counts_raise` FAIL.

- [ ] **Step 3: Replace `src/zfc/zfoldchange.py` entirely**

Write the complete new file (remove `from .statsfunc import ecdf`; keep everything else):

```python
####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn import linear_model
from .statsfunc import df_normalization
from .statsfunc import df_smallcount
from .statsfunc import p_adjust
from .statsfunc import df_robust_rank_aggregation
from .statsfunc import df_mean_rank_aggregation


def zfoldchange(data,
                top_n_sgrna=None,
                top_n_gene=None,
                iteration=100,
                normalization='total',
                seed=0):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in data.columns, 'data should have column: {}'.format(x)

    data = data.copy()

    # Validate identifiers and counts.
    for col in ['gene', 'guide', 'barcode']:
        if data[col].isna().any():
            raise ValueError(
                "column '{}' contains missing identifiers".format(col)
            )
    for col in ['ctrl', 'exp']:
        values = pd.to_numeric(data[col], errors='raise')
        if (values < 0).any():
            raise ValueError(
                "column '{}' contains negative counts".format(col)
            )
        if (values > 0).sum() == 0:
            raise ValueError(
                "column '{}' has no positive counts".format(col)
            )
        data[col] = values.fillna(0.0)

    # ------------------
    # Step 1: Normalization of raw counts
    norm = df_normalization(
        data[['ctrl', 'exp']],
        normalization
    )
    for a in [0.05, 0.1, 0.15]:
        smallcount = df_smallcount(norm, a)
        if len(smallcount) > 0:
            break

    if len(smallcount) == 0:
        smallcount = np.array(
            [data['ctrl'][data['ctrl'] > 0].quantile(0.15)]
        )

    norm = norm + smallcount.mean()
    if (norm['ctrl'] <= 0).any():
        raise ValueError(
            'ctrl counts are not positive after normalization'
        )
    bar_df = pd.concat(
        [
            data[['gene', 'guide', 'barcode']],
            norm
        ],
        axis=1, sort=False
    )
    del norm

    # ------------------
    # Step 2: Calculate fold change
    bar_df['fc'] = bar_df['exp'] / bar_df['ctrl']
    bar_df['lfc'] = np.log2(bar_df['fc'])

    # ------------------
    # Step 3: Calculate fold change std
    lfc_q1 = bar_df['lfc'].quantile(0.05)
    lfc_q2 = bar_df['lfc'].quantile(0.95)
    ctrl_q1 = bar_df['ctrl'].quantile(0.025)
    ctrl_q2 = bar_df['ctrl'].quantile(0.95)

    model_data = bar_df.loc[
        (bar_df['lfc'] >= lfc_q1) & (bar_df['lfc'] <= lfc_q2) &
        (bar_df['ctrl'] >= ctrl_q1) & (bar_df['ctrl'] <= ctrl_q2)
    ]

    bins = pd.cut(model_data['ctrl'], 200).astype(str)

    train_data = pd.DataFrame(
        {
            'lfcstd': model_data.groupby(bins)['lfc'].std(),
            'ctrlmean': model_data.groupby(bins)['ctrl'].mean()
        }
    )

    train_data = train_data.loc[~train_data['lfcstd'].isnull()]

    reg = linear_model.LinearRegression()
    reg.fit(
        np.expand_dims(train_data['ctrlmean'].values, 1),
        train_data['lfcstd'].values
    )

    bar_df = bar_df.assign(
        lfc_std=bar_df['ctrl'] * reg.coef_[0] + reg.intercept_
    )

    bar_df.loc[bar_df['lfc_std'] <= 0, 'lfc_std'] = bar_df.loc[
        bar_df['lfc_std'] > 0,
        'lfc_std'
    ].median()
    bar_df.loc[
        bar_df['ctrl'] > train_data['ctrlmean'].max(),
        'lfc_std'
    ] = bar_df['lfc_std'].median()

    # ------------------
    # Step 4: Calculate raw zlfc using lfc and lfc_std
    bar_df.loc[:, 'zlfc'] = bar_df['lfc'] / bar_df['lfc_std']

    bar_df['rank_down'] = bar_df['zlfc'].rank(
        method='average', ascending=True
    ) / len(bar_df['zlfc'])

    bar_df['rank_up'] = bar_df['zlfc'].rank(
        method='average', ascending=False
    ) / len(bar_df['zlfc'])

    bar_df.set_index(['gene', 'guide', 'barcode'], inplace=True)

    # Calculate P value using normal distribution (vectorized one-sided p).
    bar_p = stats.norm.cdf(bar_df['zlfc'])
    bar_df['p'] = np.minimum(bar_p, 1 - bar_p)

    rng = (
        np.random.default_rng()
        if (seed is None or seed < 0)
        else np.random.default_rng(seed)
    )

    # ------------------
    # Step 5: Calculate sgRNA mean zscore of fold change
    sg_groups = bar_df.groupby(level=[0, 1])['zlfc']
    sg_df = pd.DataFrame(
        {
            'zlfc': sg_groups.mean(),
            'count': sg_groups.count()
        }
    )
    sg_codes = sg_groups.ngroup().to_numpy()
    sg_counts = np.bincount(sg_codes)
    sg_uniq = np.unique(sg_counts)
    sg_nulls = {int(c): [] for c in sg_uniq}
    sg_vals = bar_df['zlfc'].to_numpy()

    for _ in range(iteration):
        perm = rng.permutation(sg_vals)
        means = np.bincount(
            sg_codes, weights=perm, minlength=len(sg_df)
        ) / sg_counts
        for c in sg_uniq:
            sg_nulls[int(c)].append(means[sg_counts == c])
    for c in sg_uniq:
        sg_nulls[int(c)] = np.sort(np.concatenate(sg_nulls[int(c)]))

    sg_z = sg_df['zlfc'].to_numpy()
    sg_p = np.empty(len(sg_df), dtype=float)
    for c in sg_uniq:
        mask = sg_counts == c
        null = sg_nulls[int(c)]
        ecdf_vals = np.searchsorted(null, sg_z[mask], side='left') / null.size
        sg_p[mask] = np.minimum(ecdf_vals, 1 - ecdf_vals)
    sg_df['p'] = sg_p
    sg_df['p_adj'] = p_adjust(sg_df['p'], 'BH')
    del sg_nulls

    # -------------------
    # Step 6: Calculate Gene mean zscore of fold change
    g_groups = sg_df.groupby(level=0)['zlfc']
    g_df = pd.DataFrame(
        {
            'zlfc': g_groups.mean(),
            'count': g_groups.count()
        }
    )
    g_codes = g_groups.ngroup().to_numpy()
    g_counts = np.bincount(g_codes)
    g_uniq = np.unique(g_counts)
    g_nulls = {int(c): [] for c in g_uniq}
    g_vals = sg_df['zlfc'].to_numpy()

    for _ in range(iteration):
        perm = rng.permutation(g_vals)
        means = np.bincount(
            g_codes, weights=perm, minlength=len(g_df)
        ) / g_counts
        for c in g_uniq:
            g_nulls[int(c)].append(means[g_counts == c])
    for c in g_uniq:
        g_nulls[int(c)] = np.sort(np.concatenate(g_nulls[int(c)]))

    g_z = g_df['zlfc'].to_numpy()
    g_p = np.empty(len(g_df), dtype=float)
    for c in g_uniq:
        mask = g_counts == c
        null = g_nulls[int(c)]
        ecdf_vals = np.searchsorted(null, g_z[mask], side='left') / null.size
        g_p[mask] = np.minimum(ecdf_vals, 1 - ecdf_vals)
    g_df['p'] = g_p
    g_df['p_adj'] = p_adjust(g_df['p'], 'BH')
    del g_nulls

    # ------------------
    # Step 7: Rank aggregation

    # sgRNA Rank
    if top_n_sgrna is None:
        top_n_sgrna = int(
            bar_df.groupby(level=[0, 1])['zlfc'].count().median()
        )
    sg_b_down = bar_df[['rank_down']].copy()
    sg_b_down.loc[:, 'groupid'] = sg_b_down.groupby(
        level=[0, 1]
    )['rank_down'].rank(
        method='first', ascending=True
    ).astype(int)
    sg_b_down.reset_index(drop=False, inplace=True)
    del sg_b_down['barcode']
    sg_b_down.set_index(['gene', 'guide', 'groupid'], inplace=True)
    sg_b_down = sg_b_down.unstack(level=2)
    sg_b_down.columns = sg_b_down.columns.levels[1]
    sg_b_down = sg_b_down[list(range(1, top_n_sgrna + 1))]

    sg_b_up = bar_df[['rank_up']].copy()
    sg_b_up.loc[:, 'groupid'] = sg_b_up.groupby(
        level=[0, 1]
    )['rank_up'].rank(
        method='first', ascending=True
    ).astype(int)
    sg_b_up.reset_index(drop=False, inplace=True)
    del sg_b_up['barcode']
    sg_b_up.set_index(['gene', 'guide', 'groupid'], inplace=True)
    sg_b_up = sg_b_up.unstack(level=2)
    sg_b_up.columns = sg_b_up.columns.levels[1]
    sg_b_up = sg_b_up[list(range(1, top_n_sgrna + 1))]

    sg_df.loc[:, 'RRA_Score_down'] = df_robust_rank_aggregation(sg_b_down)
    sg_df.loc[:, 'RRA_Score_down_adj'] = p_adjust(
        sg_df['RRA_Score_down'], 'BH'
    )
    sg_df.loc[:, 'RRA_Score_up'] = df_robust_rank_aggregation(sg_b_up)
    sg_df.loc[:, 'RRA_Score_up_adj'] = p_adjust(
        sg_df['RRA_Score_up'], 'BH'
    )

    sg_df.loc[:, 'Mean_Rank_down'] = df_mean_rank_aggregation(sg_b_down)
    sg_df.loc[:, 'Mean_Rank_up'] = df_mean_rank_aggregation(sg_b_up)

    # gene Rank
    if top_n_gene is None:
        top_n_gene = int(
            bar_df.groupby(level=0)['zlfc'].count().median()
        )
    g_b_down = bar_df[['rank_down']].copy()
    g_b_down.loc[:, 'groupid'] = g_b_down.groupby(
        level=0
    )['rank_down'].rank(
        method='first', ascending=True
    ).astype(int)
    g_b_down.reset_index(drop=False, inplace=True)
    del g_b_down['barcode']
    del g_b_down['guide']
    g_b_down.set_index(['gene', 'groupid'], inplace=True)
    g_b_down = g_b_down.unstack(level=1)
    g_b_down.columns = g_b_down.columns.levels[1]
    g_b_down = g_b_down[list(range(1, top_n_gene + 1))]

    g_b_up = bar_df[['rank_up']].copy()
    g_b_up.loc[:, 'groupid'] = g_b_up.groupby(
        level=0
    )['rank_up'].rank(
        method='first', ascending=True
    ).astype(int)
    g_b_up.reset_index(drop=False, inplace=True)
    del g_b_up['barcode']
    del g_b_up['guide']
    g_b_up.set_index(['gene', 'groupid'], inplace=True)
    g_b_up = g_b_up.unstack(level=1)
    g_b_up.columns = g_b_up.columns.levels[1]
    g_b_up = g_b_up[list(range(1, top_n_gene + 1))]

    g_df.loc[:, 'RRA_Score_down'] = df_robust_rank_aggregation(g_b_down)
    g_df.loc[:, 'RRA_Score_down_adj'] = p_adjust(
        g_df['RRA_Score_down'], 'BH'
    )
    g_df.loc[:, 'RRA_Score_up'] = df_robust_rank_aggregation(g_b_up)
    g_df.loc[:, 'RRA_Score_up_adj'] = p_adjust(
        g_df['RRA_Score_up'], 'BH'
    )

    g_df.loc[:, 'Mean_Rank_down'] = df_mean_rank_aggregation(g_b_down)
    g_df.loc[:, 'Mean_Rank_up'] = df_mean_rank_aggregation(g_b_up)

    return (bar_df, sg_df, g_df, train_data, reg)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `uv run pytest tests/test_zfoldchange.py -v`
Expected: all PASS (including the existing end-to-end test).

- [ ] **Step 5: Commit**

```bash
git add src/zfc/zfoldchange.py tests/test_zfoldchange.py
git commit -m "feat: seed, validation, and vectorized null distribution in zfoldchange"
```

---

### Task 6: `programio.py` — vectorized index + clear missing-column error

**Files:**
- Modify: `tests/test_io.py`
- Modify: `src/zfc/programio.py:16-27`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_io.py`:

```python
def test_read_raw_count_missing_column_raises(tmp_path):
    path = tmp_path / 'bad.tsv'
    path.write_text('gene\tguide\tctrl\texp\nA1BG\tG1\t100\t50\n')
    with pytest.raises(ValueError, match='barcode'):
        read_raw_count(str(path))
```

(Add `import pytest` at the top of `tests/test_io.py`.)

- [ ] **Step 2: Run test to verify it fails**

Run: `uv run pytest tests/test_io.py::test_read_raw_count_missing_column_raises -v`
Expected: FAIL (current code raises `KeyError`, not `ValueError`).

- [ ] **Step 3: Implement validation + vectorized index**

In `src/zfc/programio.py`, replace the end of `read_raw_count` (lines 16–27) with:

```python
    required = [
        gene_colname, guide_colname, barcode_colname,
        ctrl_colname, exp_colname,
    ]
    missing = [c for c in required if c not in data.columns]
    if missing:
        raise ValueError(
            'input is missing required column(s): {}'.format(', '.join(missing))
        )
    data = data[required]
    data.columns = ['gene', 'guide', 'barcode', 'ctrl', 'exp']
    data.index = (
        data['gene'] + '_' + data['guide'] + '_' + data['barcode']
    )
    return data
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_io.py -v`
Expected: PASS (existing `test_read_raw_count` index assertion `['A1BG_G1_B1', ...]` still holds — the vectorized join is string-identical).

- [ ] **Step 5: Commit**

```bash
git add src/zfc/programio.py tests/test_io.py
git commit -m "perf: vectorize index join and clarify missing-column error"
```

---

### Task 7: `cli.py` — `--seed` and input-existence error

**Files:**
- Modify: `tests/test_cli.py`
- Modify: `src/zfc/cli.py`

- [ ] **Step 1: Write the failing tests**

Replace `tests/test_cli.py` with:

```python
import numpy as np
import pytest

from zfc.cli import main


def test_help_exits_zero():
    with pytest.raises(SystemExit) as exc_info:
        main(['--help'])
    assert exc_info.value.code == 0


def test_version_exits_zero():
    with pytest.raises(SystemExit) as exc_info:
        main(['--version'])
    assert exc_info.value.code == 0


def _write_counts(path):
    rng = np.random.default_rng(0)
    lines = ['gene\tguide\tbarcode\tctrl\texp']
    for g in range(50):
        for s in range(8):
            for b in range(4):
                ctrl = float(rng.integers(50, 4000))
                exp = ctrl * float(rng.uniform(0.3, 2.0))
                lines.append(
                    'gene{0}\tguide{0}_{1}\tbar{0}_{1}_{2}\t{3}\t{4}'.format(
                        g, s, b, ctrl, exp
                    )
                )
    path.write_text('\n'.join(lines) + '\n')


def test_seed_flag_produces_reproducible_output(tmp_path):
    counts = tmp_path / 'counts.tsv'
    _write_counts(counts)
    main(['-i', str(counts), '--seed', '3', '--outprefix', str(tmp_path / 'a')])
    first = (tmp_path / 'a_sgrna.txt').read_text()

    main(['-i', str(counts), '--seed', '3', '--outprefix', str(tmp_path / 'b')])
    second = (tmp_path / 'b_sgrna.txt').read_text()

    assert first == second


def test_missing_input_errors(tmp_path):
    with pytest.raises(SystemExit) as exc_info:
        main(['-i', str(tmp_path / 'nope.tsv'), '--outprefix', str(tmp_path / 'o')])
    assert exc_info.value.code == 2
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `uv run pytest tests/test_cli.py -v`
Expected: `test_seed_flag_produces_reproducible_output` FAIL (argparse rejects `--seed` with exit code 2); `test_missing_input_errors` may PASS or raise FileNotFoundError (either way the `--seed` test gates the task).

- [ ] **Step 3: Edit `src/zfc/cli.py`**

Add `import os` after `import argparse`.

Add the `--seed` argument (after the `--null-iteration` block):

```python
    parser.add_argument(
        '--seed',
        type=int,
        default=0,
        help='Random seed for the null-distribution permutations. ' +
        'Default 0 (fixed) for reproducible output; ' +
        'use -1 for nondeterministic output.'
    )
```

Add the input check after `args = vars(parser.parse_args(argv))`:

```python
    input_path = args['input']
    if not os.path.isfile(input_path):
        parser.error('input file not found: {}'.format(input_path))
```

Change `read_raw_count(args['input'], ...)` to `read_raw_count(input_path, ...)`.

Add `seed=args['seed']` to the `zfc.zfoldchange(...)` call.

- [ ] **Step 4: Run tests to verify they pass**

Run: `uv run pytest tests/test_cli.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/zfc/cli.py tests/test_cli.py
git commit -m "feat: add --seed and clearer input-file error to CLI"
```

---

### Task 8: `plot.py` — scalar colors, close figures, no caller mutation

**Files:**
- Create: `tests/test_plot.py`
- Modify: `src/zfc/plot.py`

- [ ] **Step 1: Write the smoke test**

Create `tests/test_plot.py`:

```python
import numpy as np
import pandas as pd

from zfc import zfoldchange, zfc_plot


def test_zfc_plot_writes_outputs(tmp_path):
    rng = np.random.default_rng(3)
    rows = []
    for g in range(40):
        for s in range(6):
            for b in range(4):
                ctrl = float(rng.integers(50, 4000))
                exp = ctrl * float(rng.uniform(0.3, 2.0))
                rows.append((f'gene{g}', f'guide{g}_{s}', f'bar{g}_{s}_{b}', ctrl, exp))
    data = pd.DataFrame(rows, columns=['gene', 'guide', 'barcode', 'ctrl', 'exp'])

    bar, sg, g, train, reg = zfoldchange(data, iteration=5, seed=0)
    out = str(tmp_path / 'out')
    zfc_plot(out, data, bar, sg, g, train, reg)

    assert (tmp_path / 'out_lorenz.png').exists()
    assert (tmp_path / 'out_gene_zlfc_hist.png').exists()
    assert (tmp_path / 'out_sgrna_zlfc_rra_fdr_scatter.pdf').exists()
```

- [ ] **Step 2: Run test to verify it passes on current code**

Run: `uv run pytest tests/test_plot.py -v`
Expected: PASS (current plotting works; this test guards the refactor against regressions).

- [ ] **Step 3: Apply the plot fixes**

In `src/zfc/plot.py`, add the helper after the `gini` import:

```python
def _save(fig, outprefix, name):
    fig.savefig('_'.join([outprefix, name + '.png']))
    fig.savefig('_'.join([outprefix, name + '.pdf']))
    plt.close(fig)
```

Replace every `c=['black'] * len(...), alpha=0.1, edgecolor=['none'] * len(...)` pattern with `c='black', alpha=0.1, edgecolor='none'` (7 sites: `normcount_scatter`, `lfc_normcount_scatter` ×2, `zlfc_normcount_scatter` ×2, `lfcstd_lfc_scatter`, `zlfc_rra_scatter`, `zlfc_rrafdr_scatter`).

In `zlfc_rra_scatter`, replace the mutating lines:

```python
    data.loc[:, 'RRA'] = data[['RRA_Score_down', 'RRA_Score_up']].min(axis=1)
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        data['zlfc'], np.log10(data['RRA']) * -1,
        c=['black'] * len(data['zlfc']), alpha=0.1,
        edgecolor=['none'] * len(data['zlfc'])
    )
```
with:
```python
    rra = data[['RRA_Score_down', 'RRA_Score_up']].min(axis=1)
    fig, axes = plt.subplots(1, 1)
    axes.scatter(data['zlfc'], np.log10(rra) * -1, c='black', alpha=0.1, edgecolor='none')
```

Do the same in `zlfc_rrafdr_scatter` (use `RRA_Score_down_adj`/`RRA_Score_up_adj`), and rename its local to `rra`.

Replace the body of `zfc_plot` (from `# Lorenz Curve` through the last `savefig`) with:

```python
def zfc_plot(outprefix, inputdata,
             barresult, sgresult, gresult,
             train_data, lm):
    _save(lorenz_curve(inputdata), outprefix, 'lorenz')
    _save(counts_boxplot(inputdata, barresult), outprefix, 'counts_boxplot')
    _save(normcount_scatter(barresult), outprefix, 'normcount_scatter')
    _save(lfc_normcount_scatter(barresult), outprefix, 'lfc_normcount_scatter')
    _save(zlfc_normcount_scatter(barresult), outprefix, 'zlfc_normcount_scatter')
    _save(lfcstd_lfc_scatter(barresult), outprefix, 'lfc_std_scatter')
    _save(zlfc_hist(sgresult), outprefix, 'sgrna_zlfc_hist')
    _save(zlfc_hist(gresult), outprefix, 'gene_zlfc_hist')
    _save(zlfc_rra_scatter(sgresult), outprefix, 'sgrna_zlfc_rra_scatter')
    _save(zlfc_rra_scatter(gresult), outprefix, 'gene_zlfc_rra_scatter')
    _save(zlfc_rrafdr_scatter(sgresult), outprefix, 'sgrna_zlfc_rra_fdr_scatter')
    _save(zlfc_rrafdr_scatter(gresult), outprefix, 'gene_zlfc_rra_fdr_scatter')
    _save(lfcstd_ctrlmean_figure(train_data, lm), outprefix, 'lfcstd_ctrlmean')
```

(Optional clarity: rename the second parameter of `counts_boxplot` from `sgresult` to `barresult` — call sites are positional, so no behavior change.)

- [ ] **Step 4: Run the test to verify it still passes**

Run: `uv run pytest tests/test_plot.py -v`
Expected: PASS, and all 13 × 2 output files are written.

- [ ] **Step 5: Commit**

```bash
git add src/zfc/plot.py tests/test_plot.py
git commit -m "perf: scalar plot colors and close figures; stop mutating results in plots"
```

---

### Task 9: Bump version and update docs

**Files:**
- Modify: `src/zfc/__init__.py:17`
- Modify: `README.md`
- Modify: `docs/algorithm.md`

- [ ] **Step 1: Bump the version**

In `src/zfc/__init__.py`, change `__version__ = '0.1.7'` to `__version__ = '0.2.0'`.

- [ ] **Step 2: Update README normalization + seed help**

In `README.md`, change the `--normalization` help block (lines 68–73) to:

```text
  --normalization  {total,median,upper_quartile,median_ratio,none}
                        Normalization of raw count data, default is total.
                        Support method: total (Total sum normalization);
                        median (Median normalization); upper_quartile (Upper
                        quartile normalization (0.75)); median_ratio (Median
                        ratio normalization); none (No normalization).
  --seed SEED            Random seed for the null-distribution permutations.
                        Default 0 (fixed) for reproducible output; use -1 for
                        nondeterministic output.
```

(Also add the same `--seed` line in the `--null-iteration` vicinity of the help snapshot so the sample output matches the new CLI. Verify with `zfc --help` and reconcile any wording.)

- [ ] **Step 3: Update `docs/algorithm.md` "Known issues"**

In `docs/algorithm.md`, replace the `### Known bug: upper_quartile branch never runs` section and the `## Known issues` list (items 1–5) to reflect the resolution:
- upper_quartile now routes to `df_quantile_normfactor`;
- `median_ratio` warns + falls back to `total` (no `exit()`);
- `df_rank_order_aggregation` deleted (NumPy 2 `np.math`).
- reproducibility: `--seed` added (default 0; `-1` = nondeterministic);
- div-by-zero: input validation + post-normalization positive-ctrl guard.

Mark each as fixed with a one-line note (keep the doc's "fix deliberately" tone).

- [ ] **Step 4: Lint and test**

Run: `uvx ruff check` then `uv run pytest`
Expected: clean and green.

- [ ] **Step 5: Commit**

```bash
git add src/zfc/__init__.py README.md docs/algorithm.md
git commit -m "docs: bump to 0.2.0 and resolve known-issue notes"
```

---

### Task 10: Final verification and benchmark

**Files:** none (verification only)

- [ ] **Step 1: Full lint + test**

Run: `uvx ruff check && uv run pytest`
Expected: lint clean, all tests green.

- [ ] **Step 2: Benchmarks**

Run: `.venv/bin/python -c "import time; from zfc import read_raw_count, zfoldchange; d = read_raw_count('tmp/REP_A549_High.txt'); t = time.time(); zfoldchange(d, iteration=3, seed=0); print('iter=3 %.2fs' % (time.time() - t))"`

Expected: **~2–3s** (baseline was 19.8s at iteration=3).

Run: `.venv/bin/python -c "import time; from zfc import read_raw_count, zfoldchange; d = read_raw_count('tmp/REP_A549_High.txt'); t = time.time(); zfoldchange(d, iteration=100, seed=0); print('iter=100 %.2fs' % (time.time() - t))"`

Expected: **~4–6s** (measured baseline 82.2s at iteration=100; the remaining time is dominated by the rank-matrix build, not the null loop).

- [ ] **Step 3: Determinism + contract spot-check (one file, real data)**

Run: `.venv/bin/python -c "from zfc import read_raw_count, zfoldchange; d = read_raw_count('tmp/REP_A549_High.txt'); a = zfoldchange(d, iteration=5, seed=0)[1]; b = zfoldchange(d, iteration=5, seed=0)[1]; print('sgrna p identical:', (a['p'] == b['p']).all()); print('sgrna columns:', list(a.columns))"`

Expected: `sgrna p identical: True`; columns match the contract (`zlfc, count, p, p_adj, RRA_Score_down, RRA_Score_down_adj, RRA_Score_up, RRA_Score_up_adj, Mean_Rank_down, Mean_Rank_up`).

- [ ] **Step 4: Record results**

Record the before/after timings and the determinism confirmation in the commit message or PR description; note the ≥5× speedup.

- [ ] **Step 5: Final commit (if any benchmarking notes added)**

```bash
git add -A
git commit -m "bench: record optimization timings"   # only if notes were added
```

---

## Self-Review

**Spec coverage** (from `docs/specs/2026-08-21-zfc-optimization-design.md`):
- §1 RRA vectorization → Task 1; §2 upper_quartile/median_ratio/delete → Tasks 2–4; §3 seed/null/guards/validation → Task 5; §4 programio → Task 6; §5 cli → Task 7; §6 plot → Task 8; §7 version → Task 9; tests → tests checked in each task; success criteria/before-after → Task 10. ✅ No gap.

**Placeholder scan:** no TBD/TODO; every code step shows full code and every command has an expected result. ✅

**Type/signature consistency:** `zfoldchange(..., seed=0)` (Task 5) matches `--seed` default `0` (Task 7) and `seed=args['seed']`; `df_robust_rank_aggregation` returns `pd.Series` (Task 1) consumed as before in Task 5; `np.minimum` p-reflection consistent across `bar_df` and the null steps. ✅