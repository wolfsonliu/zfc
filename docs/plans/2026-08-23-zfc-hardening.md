# ZFC Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Harden ZFC with robustness fixes, an output-preserving refactor of `zfoldchange`, type annotations, plot cleanup, CI, and README fixes — without changing any numeric result.

**Architecture:** Strictly output-preserving. A golden-parity test (`tests/test_parity.py`) locks the current `0.2.1` barcode/sgRNA/gene tables as SHA-256 digests before any refactor; every later change must keep that test green. Existing tests plus new behavior tests guard each step. No performance work, no statistical-convention changes, no new features.

**Tech Stack:** Python >= 3.11, pandas/numpy/scipy/scikit-learn/matplotlib, uv + ruff + pytest (existing), GitHub Actions (new).

**Spec:** `docs/specs/2026-08-23-zfc-hardening-design.md`.

---

## Environment / setup note (read first)

The project's canonical commands are `uv run pytest` and `uvx ruff check`
(`docs/development.md`). If either fails with `Could not acquire lock` /
`Read-only file system ... /home/zhiheng/.cache/uv`, the uv cache dir is not
writable in this workspace. Prefix those commands with a writable cache:

```bash
export UV_CACHE_DIR=/tmp/uv-cache
```

All task commands below omit the prefix for readability; add it only if you
hit the lock error.

---

## Task 1: Golden-parity regression test (locks 0.2.1 output)

This is a *characterization* test, not a red-green test: it asserts current
behavior and is expected to PASS immediately. It becomes the guard for Tasks
2–6.

**Files:**
- Create: `tests/test_parity.py`

- [ ] **Step 1: Write the parity test**

```python
import hashlib

import numpy as np
import pandas as pd

from zfc.zfoldchange import zfoldchange


# Captured from zfc 0.2.1 running on _make_input() with iteration=10, seed=0.
# If these ever change, STOP: the change is a numeric-output regression.
BAR_SHA256 = '90204f84d1d659048837cadec2ec49af8cec2e1637334b25a7fbe8b1352c0075'
SG_SHA256 = '4cfc08d47990992be963f3f8d57bb4bb44535fb8dff821bdace0f296880eb013'
G_SHA256 = '21848c5341e2da402f278e7f9b140cdac843dd600003acaefebf7eae89d726c8'


def _make_input(n_genes=20, n_guides=6, n_barcodes=3, seed=7):
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        for s in range(n_guides):
            for b in range(n_barcodes):
                ctrl = float(rng.integers(50, 4000))
                exp = ctrl * float(rng.uniform(0.3, 2.0))
                rows.append((
                    'gene{}'.format(g), 'guide{}_{}'.format(g, s),
                    'bar{}_{}_{}'.format(g, s, b), ctrl, exp,
                ))
    return pd.DataFrame(
        rows, columns=['gene', 'guide', 'barcode', 'ctrl', 'exp']
    )


def _digest(frame):
    payload = frame.to_csv(float_format='%.17g', sep='\t').encode('utf-8')
    return hashlib.sha256(payload).hexdigest()


def test_golden_output_parity():
    bar, sg, g, _, _ = zfoldchange(_make_input(), iteration=10, seed=0)
    assert _digest(bar) == BAR_SHA256
    assert _digest(sg) == SG_SHA256
    assert _digest(g) == G_SHA256
```

- [ ] **Step 2: Run it — must pass**

Run: `uv run pytest tests/test_parity.py -v`
Expected: `1 passed` (this locks the baseline; it is NOT supposed to fail)

- [ ] **Step 3: Commit**

```bash
git add tests/test_parity.py
git commit -m "test: add golden-parity regression test for result tables"
```

---

## Task 2: Make `--input` required (no more TypeError on bare `zfc`)

**Files:**
- Modify: `src/zfc/cli.py:14-23` (the `-i/--input` argument)
- Test: `tests/test_cli.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_cli.py`:

```python
def test_missing_input_flag_errors():
    with pytest.raises(SystemExit) as exc_info:
        main([])  # no -i at all
    assert exc_info.value.code == 2
```

- [ ] **Step 2: Run it — must FAIL (proves the bug)**

Run: `uv run pytest tests/test_cli.py::test_missing_input_flag_errors -v`
Expected: FAIL — currently `os.path.isfile(None)` raises `TypeError`, so
`pytest.raises(SystemExit)` does not fire.

- [ ] **Step 3: Fix the argument**

In `src/zfc/cli.py`, change:

```python
    parser.add_argument(
        '-i', '--input',
        action='store',
        help='Raw count table with header should include' +
```

to:

```python
    parser.add_argument(
        '-i', '--input',
        action='store',
        required=True,
        help='Raw count table with header should include' +
```

(`required=True` makes argparse exit code 2 with a clear message before any
`os.path.isfile` call.)

- [ ] **Step 4: Run the CLI tests — must pass**

Run: `uv run pytest tests/test_cli.py -v`
Expected: all pass (new test + existing `test_missing_input_errors`,
`test_help_exits_zero`, `test_version_exits_zero`, `test_seed_flag...`)

- [ ] **Step 5: Commit**

```bash
git add src/zfc/cli.py tests/test_cli.py
git commit -m "fix: require --input in the CLI"
```

---

## Task 3: `zfoldchange` missing-column `assert` → `ValueError`

**Files:**
- Modify: `src/zfc/zfoldchange.py:26-27`
- Test: `tests/test_zfoldchange.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_zfoldchange.py`:

```python
def test_missing_column_raises_value_error():
    data = _make_input()
    data = data.drop(columns=['barcode'])
    with pytest.raises(ValueError, match='barcode'):
        zfoldchange(data, iteration=5, seed=0)
```

- [ ] **Step 2: Run it — must FAIL**

Run: `uv run pytest tests/test_zfoldchange.py::test_missing_column_raises_value_error -v`
Expected: FAIL — current code raises `AssertionError`, not `ValueError`.

- [ ] **Step 3: Replace the assert**

In `src/zfc/zfoldchange.py`, replace lines 26-27:

```python
    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in data.columns, 'data should have column: {}'.format(x)
```

with:

```python
    missing = [
        c for c in ['gene', 'guide', 'barcode', 'ctrl', 'exp']
        if c not in data.columns
    ]
    if missing:
        raise ValueError(
            'data should have column(s): {}'.format(', '.join(missing))
        )
```

- [ ] **Step 4: Run the zfoldchange tests — must pass**

Run: `uv run pytest tests/test_zfoldchange.py -v`
Expected: all pass (new test + existing end-to-end, determinism, robustness tests).

- [ ] **Step 5: Commit**

```bash
git add src/zfc/zfoldchange.py tests/test_zfoldchange.py
git commit -m "fix: raise ValueError for missing zfoldchange columns"
```

---

## Task 4: Split `zfoldchange` into step helpers (output-preserving)

This replaces `zfoldchange.py` with an equivalently-structured module. Every
numeric operation is moved verbatim (same ordering, same RNG draws), so Tasks
1's parity digest must still pass.

**Files:**
- Modify: `src/zfc/zfoldchange.py` (whole file)

- [ ] **Step 1: Replace the file**

Write the full file (this is the entire target content):

```python
####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.stats as stats
from sklearn import linear_model

from .statsfunc import df_normalization
from .statsfunc import df_smallcount
from .statsfunc import p_adjust
from .statsfunc import df_robust_rank_aggregation
from .statsfunc import df_mean_rank_aggregation

REQUIRED_COLUMNS = ['gene', 'guide', 'barcode', 'ctrl', 'exp']


def _validate_input(data: pd.DataFrame) -> pd.DataFrame:
    """Validate a raw count frame and return a cleaned copy."""
    missing = [c for c in REQUIRED_COLUMNS if c not in data.columns]
    if missing:
        raise ValueError(
            'data should have column(s): {}'.format(', '.join(missing))
        )

    data = data.copy()

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

    return data


def _normalize(data: pd.DataFrame, normalization: str) -> pd.DataFrame:
    """Apply normalization + small-count bump; return the barcode frame."""
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
    return pd.concat(
        [
            data[['gene', 'guide', 'barcode']],
            norm
        ],
        axis=1, sort=False
    )


def _fold_change(bar_df: pd.DataFrame) -> pd.DataFrame:
    """Add fc and lfc columns."""
    bar_df['fc'] = bar_df['exp'] / bar_df['ctrl']
    bar_df['lfc'] = np.log2(bar_df['fc'])
    return bar_df


def _lfc_std_model(bar_df: pd.DataFrame):
    """Fit the lfc-std linear model and add lfc_std (returns bar_df, train, reg)."""
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

    return bar_df, train_data, reg


def _barcode_zlfc_rank(bar_df: pd.DataFrame) -> pd.DataFrame:
    """Add zlfc, rank_down/up, and the barcode-level normal p."""
    bar_df.loc[:, 'zlfc'] = bar_df['lfc'] / bar_df['lfc_std']

    bar_df['rank_down'] = bar_df['zlfc'].rank(
        method='average', ascending=True
    ) / len(bar_df['zlfc'])

    bar_df['rank_up'] = bar_df['zlfc'].rank(
        method='average', ascending=False
    ) / len(bar_df['zlfc'])

    bar_df.set_index(['gene', 'guide', 'barcode'], inplace=True)

    bar_p = stats.norm.cdf(bar_df['zlfc'])
    bar_df['p'] = np.minimum(bar_p, 1 - bar_p)
    return bar_df


def _empirical_p(groups, values: pd.Series, iteration: int, rng) -> pd.DataFrame:
    """Per-group empirical p-values via permutation (sgRNA and gene levels)."""
    group_df = pd.DataFrame(
        {
            'zlfc': groups.mean(),
            'count': groups.count()
        }
    )
    codes = groups.ngroup().to_numpy()
    counts = np.bincount(codes)
    uniq = np.unique(counts)
    nulls = {int(c): [] for c in uniq}
    vals = values.to_numpy()

    for _ in range(iteration):
        perm = rng.permutation(vals)
        means = np.bincount(
            codes, weights=perm, minlength=len(group_df)
        ) / counts
        for c in uniq:
            nulls[int(c)].append(means[counts == c])
    for c in uniq:
        nulls[int(c)] = np.sort(np.concatenate(nulls[int(c)]))

    z = group_df['zlfc'].to_numpy()
    p = np.empty(len(group_df), dtype=float)
    for c in uniq:
        mask = counts == c
        null = nulls[int(c)]
        ecdf_vals = np.searchsorted(null, z[mask], side='left') / null.size
        p[mask] = np.minimum(ecdf_vals, 1 - ecdf_vals)

    group_df['p'] = p
    group_df['p_adj'] = p_adjust(group_df['p'], 'BH')
    return group_df


def _rank_aggregation(bar_df: pd.DataFrame, level: str, top_n: int | None) -> dict:
    """RRA and mean-rank scores for one level ('sg' or 'gene')."""
    if level == 'sg':
        group_level = [0, 1]
        index_names = ['gene', 'guide', 'groupid']
        drop = ['barcode']
        unstack_level = 2
    else:
        group_level = [0]
        index_names = ['gene', 'groupid']
        drop = ['barcode', 'guide']
        unstack_level = 1

    if top_n is None:
        top_n = int(bar_df.groupby(level=group_level)['zlfc'].count().median())

    matrices = {}
    for rank_col, ascending in [('rank_down', True), ('rank_up', False)]:
        m = bar_df[[rank_col]].copy()
        m.loc[:, 'groupid'] = m.groupby(level=group_level)[rank_col].rank(
            method='first', ascending=ascending
        ).astype(int)
        m.reset_index(drop=False, inplace=True)
        for col in drop:
            del m[col]
        m.set_index(index_names, inplace=True)
        m = m.unstack(level=unstack_level)
        m.columns = m.columns.levels[1]
        m = m[list(range(1, top_n + 1))]
        matrices[rank_col] = m

    result = {}
    result['RRA_Score_down'] = df_robust_rank_aggregation(matrices['rank_down'])
    result['RRA_Score_down_adj'] = p_adjust(result['RRA_Score_down'], 'BH')
    result['RRA_Score_up'] = df_robust_rank_aggregation(matrices['rank_up'])
    result['RRA_Score_up_adj'] = p_adjust(result['RRA_Score_up'], 'BH')
    result['Mean_Rank_down'] = df_mean_rank_aggregation(matrices['rank_down'])
    result['Mean_Rank_up'] = df_mean_rank_aggregation(matrices['rank_up'])
    return result


def zfoldchange(data: pd.DataFrame,
                top_n_sgrna: int | None = None,
                top_n_gene: int | None = None,
                iteration: int = 100,
                normalization: str = 'total',
                seed: int = 0):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    data = _validate_input(data)

    # Step 1: Normalization of raw counts
    bar_df = _normalize(data, normalization)

    # Step 2: Calculate fold change
    bar_df = _fold_change(bar_df)

    # Step 3: Calculate fold change std
    bar_df, train_data, reg = _lfc_std_model(bar_df)

    # Step 4: Calculate raw zlfc and barcode-level normal p
    bar_df = _barcode_zlfc_rank(bar_df)

    rng = (
        np.random.default_rng()
        if (seed is None or seed < 0)
        else np.random.default_rng(seed)
    )

    # Step 5: Calculate sgRNA mean zscore of fold change
    sg_df = _empirical_p(
        bar_df.groupby(level=[0, 1])['zlfc'],
        bar_df['zlfc'],
        iteration,
        rng
    )

    # Step 6: Calculate gene mean zscore of fold change
    g_df = _empirical_p(
        sg_df.groupby(level=0)['zlfc'],
        sg_df['zlfc'],
        iteration,
        rng
    )

    # Step 7: Rank aggregation
    for col, val in _rank_aggregation(bar_df, 'sg', top_n_sgrna).items():
        sg_df.loc[:, col] = val
    for col, val in _rank_aggregation(bar_df, 'gene', top_n_gene).items():
        g_df.loc[:, col] = val

    return (bar_df, sg_df, g_df, train_data, reg)
```

- [ ] **Step 2: Run parity + zfoldchange tests — must pass**

Run: `uv run pytest tests/test_parity.py tests/test_zfoldchange.py -v`
Expected: `test_golden_output_parity` PASSES (refactor is bit-for-bit identical);
all `test_zfoldchange` tests pass.

If the parity digest changes, STOP and diff against the previous version —
a helper reorders a float operation somewhere.

- [ ] **Step 3: Run the full suite + lint**

```bash
uv run pytest -q
uvx ruff check
```
Expected: all tests green; ruff reports nothing.

- [ ] **Step 4: Commit**

```bash
git add src/zfc/zfoldchange.py
git commit -m "refactor: split zfoldchange into step helpers"
```

---

## Task 5: Type-annotate `programio.py` and `statsfunc.py`

**Files:**
- Modify: `src/zfc/programio.py` (whole file)
- Modify: `src/zfc/statsfunc.py` (signatures + `from __future__ import annotations`)

- [ ] **Step 1: Replace `programio.py`**

```python
####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

from __future__ import annotations

import pandas as pd


def read_raw_count(filepath: str,
                   gene_colname: str = 'gene',
                   guide_colname: str = 'guide',
                   barcode_colname: str = 'barcode',
                   ctrl_colname: str = 'ctrl',
                   exp_colname: str = 'exp') -> pd.DataFrame:
    data = pd.read_csv(filepath, header=0, sep='\t')
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


def write_bar_result(data: pd.DataFrame, filepath: str) -> None:
    data[[
        'ctrl', 'exp', 'fc', 'lfc', 'lfc_std', 'zlfc', 'p'
    ]].sort_index(
    ).to_csv(filepath, index=True, sep='\t')


def write_sg_result(data: pd.DataFrame, filepath: str) -> None:
    data[[
        'zlfc', 'count',
        'p', 'p_adj',
        'RRA_Score_down', 'RRA_Score_down_adj',
        'RRA_Score_up', 'RRA_Score_up_adj',
        'Mean_Rank_down', 'Mean_Rank_up',
    ]].sort_values(
        'zlfc', ascending=True
    ).to_csv(filepath, index=True, sep='\t')


def write_g_result(data: pd.DataFrame, filepath: str) -> None:
    data[[
        'zlfc', 'count',
        'p', 'p_adj',
        'RRA_Score_down', 'RRA_Score_down_adj',
        'RRA_Score_up', 'RRA_Score_up_adj',
        'Mean_Rank_down', 'Mean_Rank_up',
    ]].sort_values(
        'zlfc', ascending=True
    ).to_csv(filepath, index=True, sep='\t')
```

- [ ] **Step 2: Replace `statsfunc.py` with the annotated version**

Write the full file (logic and docstrings unchanged; only the `from __future__`
import, signatures, and a citation on `bh_correction` are added).

```python
####################
# Statistics Functions
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

from __future__ import annotations

import functools
import warnings

import numpy as np
import pandas as pd
import scipy.stats as stats


def ecdf(x: np.ndarray) -> functools.partial:
    '''Calculate the empirical cumulative density function
       ====================
       Return a cumulative density function based
       on the input data.

    '''
    x = np.array(x)
    x.sort()
    n = len(x)
    y = np.concatenate([[0], np.linspace(1/n, 1, n)])

    def childfunc(sample, x, y, sorted=True):
        if not sorted:
            asort = np.argsort(x)
            x = np.take(x, asort, 0)
            y = np.take(y, asort, 0)
        idx = np.searchsorted(x, sample)
        return y[idx]

    return functools.partial(childfunc, x=x, y=y)


def bonferroni_correction(pvals: np.ndarray) -> np.ndarray:
    '''Bonferroni multitest p value correction'''

    pvals = np.asarray(pvals)
    n = len(pvals)
    pvals_corrected = pvals * float(n)
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected


def bh_correction(pvals: np.ndarray) -> np.ndarray:
    '''BH multitest p value correction.

    Benjamini, Y., and Hochberg, Y. (1995). J. Roy. Stat. Soc. B 57, 289-300.
    '''

    pvals = np.asarray(pvals)
    n = len(pvals)
    pvals_decreasing_argsort = pvals.argsort()[::-1]
    pvals_reverse_argsort = pvals_decreasing_argsort.argsort()
    pvals_corrected = np.minimum.accumulate(
        pvals[pvals_decreasing_argsort] * n / np.arange(n, 0, -1)
    )[pvals_reverse_argsort]
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected


def p_adjust(pvals, method: str = 'BH'):
    assert method in ['BH', 'Bonferroni'], 'Adjust method should be in: BH, Bonferroni'
    adjust_func = {
        'BH': bh_correction,
        'Bonferroni': bonferroni_correction
    }

    return adjust_func[method](pvals)


def df_median_ratio_normfactor(df: pd.DataFrame) -> pd.Series:
    '''
    Normalize input data indicated by the label.
    Now only median ratio normalization is available.

    Reference: Maza, E., Frasse, P., Senin, P., Bouzayen, M., and Zouine, M. (2013). Comparison of normalization methods for differential gene expression analysis in RNA-Seq experiments: A matter of relative size of studied transcriptomes. Commun Integr Biol 6, e25849-e25849.
    '''
    dfgm = np.exp(np.log(df + 1.0).sum(axis=1) / df.shape[1]) - 1.0
    dfgm[dfgm <= 0] = 1
    meanfactor = df.div(dfgm, axis=0)
    normfactor = 1 / meanfactor.median(axis=0)
    return normfactor


def df_total_count_normfactor(df: pd.DataFrame) -> pd.Series:
    colsum = df.sum(axis=0)
    normfactor = colsum.mean() / colsum
    return normfactor


def df_median_normfactor(df: pd.DataFrame) -> pd.Series:
    colmedian = df.median(axis=0)
    normfactor = colmedian.mean() / colmedian
    return normfactor


def df_quantile_normfactor(df: pd.DataFrame, quantile: float = 0.75) -> pd.Series:
    '''
    From edgeR .calcFactorQuantile
    '''
    normfactor = df_total_count_normfactor(df)
    y = df.mul(normfactor, axis=1) / df.sum(axis=0).mean()
    return y.quantile(quantile, axis=0)


def df_normalization(df: pd.DataFrame, method: str) -> pd.DataFrame:
    '''
    Normalization of raw count data.

    Normalization method:
        * "total": Total sum normalization
        * "median": Median normalization
        * "upper_quartile": Upper quartile normalization (0.75)
        * "median_ratio": Median ratio normalization
        * "none: No normalization
    '''
    assert method in [
        'total', 'median', 'upper_quartile', 'median_ratio', 'none'
    ], "Normalization method not support {0}".format(method)

    normfactor = df_total_count_normfactor(df)
    if method == 'none':
        normfactor = np.array([1]*len(normfactor))
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
    elif method == 'upper_quartile':
        normfactor = df_quantile_normfactor(df)
    elif method == 'median':
        normfactor = df_median_normfactor(df)
    elif method == 'total':
        pass
    result = df.mul(normfactor, axis=1)
    return result


def df_smallcount(df: pd.DataFrame, quantile: float = 0.1,
                  drop0: bool = True) -> np.ndarray:
    dfquantile = df.quantile(quantile, axis=0)
    smallidx = (df <= dfquantile).all(axis=1)
    smallcount = np.asarray(df[smallidx]).ravel()
    if drop0:
        smallcount = smallcount[smallcount > 0]
    return smallcount


def df_mean_rank_aggregation(rank_matrix: pd.DataFrame) -> pd.Series:
    # Kolde, R., Laur, S., Adler, P., and Vilo, J. (2012). Robust rank
    # aggregation for gene list integration and
    # meta-analysis. Bioinformatics 28, 573-580.

    rank_mean = rank_matrix.mean(axis=1, skipna=True)
    rank_num = rank_matrix.notna().sum(axis=1)
    rank_table = pd.DataFrame({'rank': rank_mean, 'n': rank_num})
    # Rank is uniform distribution, with var = 1/12(max - min)
    # Central limit theorem, var_new = var / n
    rank_table['sd'] = np.sqrt(1/12/rank_table['n'])
    rank_score = pd.Series(
        stats.norm.cdf(
            rank_table['rank'], loc=0.5, scale=rank_table['sd']
        ).clip(0, 1)
    )
    rank_score.index = rank_matrix.index
    return rank_score


def df_robust_rank_aggregation(rank_matrix: pd.DataFrame) -> pd.Series:
    # Kolde, R., Laur, S., Adler, P., and Vilo, J. (2012). Robust rank
    # aggregation for gene list integration and
    # meta-analysis. Bioinformatics 28, 573-580.
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


def gini(arr: np.ndarray) -> float:
    '''
    Using R library ineq, Gini function
    '''
    sorted_arr = arr.copy()
    sorted_arr.sort()
    n = arr.size
    G = sum(sorted_arr * (np.arange(n) + 1))
    G = 2 * G / sorted_arr.sum() - (n + 1.)
    return G/n
```

(`gini` keeps its original body here; Task 6 vectorizes it.)

- [ ] **Step 3: Run tests + lint — must pass**

```bash
uv run pytest -q
uvx ruff check
```
Expected: all green; ruff clean.

- [ ] **Step 4: Commit**

```bash
git add src/zfc/programio.py src/zfc/statsfunc.py
git commit -m "refactor: add type annotations to public API and statsfunc"
```

---

## Task 6: Plot cleanup (rename, dedupe, vectorize gini)

**Files:**
- Modify: `src/zfc/plot.py`

- [ ] **Step 1: Replace `gini`**

Replace the existing `gini` body (lines 194-203) with:

```python
def gini(arr: np.ndarray) -> float:
    '''
    Using R library ineq, Gini function
    '''
    sorted_arr = np.sort(arr)
    n = arr.size
    G = np.dot(sorted_arr, np.arange(n) + 1)
    G = 2 * G / sorted_arr.sum() - (n + 1.)
    return G / n
```

(Both the old `sort`+`sum` and this are O(n log n); `np.dot` removes the Python
loop. The value only feeds plot titles rounded to 4 decimals, never a result
table.)

- [ ] **Step 2: Rename `counts_boxplot` parameter**

Replace the function signature and its two internal references
(`sgresult['ctrl']` → `barresult['ctrl']`, `sgresult['exp']` →
`barresult['exp']`):

```python
def counts_boxplot(inputdata, barresult):
    fig, axes = plt.subplots(1, 1)
    box = axes.boxplot(
        [
            inputdata['ctrl'], inputdata['exp'],
            barresult['ctrl'], barresult['exp']
        ],
        patch_artist=True,
        sym='.'
    )
```

The `zfc_plot` call site already passes `barresult`; no change needed there.

- [ ] **Step 3: Deduplicate the 2×4 scatter builders**

Add this helper above `lfc_normcount_scatter`, and replace the two existing
function bodies with thin wrappers:

```python
def _lfc_scatter(data, value_col, ylabel):
    fig, axes = plt.subplots(
        2, 4,
        gridspec_kw={'width_ratios': [3, 1, 3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(12.8, 6.4)
    )
    axes[1, 0].scatter(data['ctrl'], data[value_col], color=BLUE, alpha=0.15, edgecolor='none')
    axes[1, 0].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel(ylabel)
    axes[0, 0].hist(data['ctrl'], bins=100, color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6)
    axes[0, 0].set_ylabel('Counts')
    axes[1, 1].hist(
        data[value_col], bins=100,
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6,
        orientation='horizontal'
    )
    axes[1, 1].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 1].set_xlabel('Counts')
    axes[0, 1].set_visible(False)

    axes[1, 2].scatter(data['exp'], data[value_col], color=BLUE, alpha=0.15, edgecolor='none')
    axes[1, 2].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 2].set_xlabel('Normalized counts of exp')
    axes[1, 2].set_ylabel(ylabel)
    axes[0, 2].hist(data['exp'], bins=100, color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6)
    axes[0, 2].set_ylabel('Counts')
    axes[1, 3].hist(
        data[value_col], bins=100,
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6,
        orientation='horizontal'
    )
    axes[1, 3].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 3].set_xlabel('Counts')
    axes[0, 3].set_visible(False)

    for ax in (axes[1, 0], axes[0, 0], axes[1, 1],
               axes[1, 2], axes[0, 2], axes[1, 3]):
        _style_ax(ax)
    return fig


def lfc_normcount_scatter(data):
    return _lfc_scatter(data, 'lfc', '$log_{2}$ fold change')


def zlfc_normcount_scatter(data):
    return _lfc_scatter(data, 'zlfc', 'Z score of $log_{2}$ fold change')
```

(`normcount_scatter` is a 2×2 layout and is left untouched.)

- [ ] **Step 4: Run tests — must pass**

Run: `uv run pytest tests/test_plot.py -v`
Expected: `test_zfc_plot_writes_outputs` passes (all figures still written).

- [ ] **Step 5: Lint + full suite**

```bash
uvx ruff check
uv run pytest -q
```
Expected: ruff clean; all green.

- [ ] **Step 6: Commit**

```bash
git add src/zfc/plot.py
git commit -m "refactor: deduplicate plot scatter builders and vectorize gini"
```

---

## Task 7: GitHub Actions CI

**Files:**
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Write the workflow**

```yaml
name: CI

on:
  push:
  pull_request:

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      fail-fast: false
      matrix:
        python-version: ["3.11", "3.12", "3.13"]
    env:
      MPLCONFIGDIR: ${{ runner.temp }}
    steps:
      - uses: actions/checkout@v4
      - name: Install uv
        uses: astral-sh/setup-uv@v5
      - name: Set up Python
        uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - name: Install dependencies
        run: uv sync
      - name: Lint
        run: uvx ruff check
      - name: Test
        run: uv run pytest
```

- [ ] **Step 2: Validate the file parses**

Run: `uvx ruff check .github/workflows/ci.yml 2>/dev/null || true`
(ruff will skip non-Python; just confirm the file is committed and YAML is
well-formed — the `on:`/`jobs:`/`steps:` keys above are the only required
structure and GitHub validates the rest on push.)

- [ ] **Step 3: Commit**

```bash
git add .github/workflows/ci.yml
git commit -m "ci: add GitHub Actions workflow"
```

> Note: you cannot run GitHub Actions locally; the workflow is exercised on the
> next `git push` / PR. The local equivalent (Tasks' `uvx ruff check` +
> `uv run pytest`) is the pre-push gate.

---

## Task 8: Sync `README.md` help and dependency sections

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Replace the Dependency section (lines 13-23)**

Replace:

```markdown
## Dependency ##

ZFC is designed with python3 and requires packages that are available
in Linux, Mac OS, and Windows.

* Python3.x
* numpy >= 1.10
* scipy >= 1.0
* pandas >= 0.16
* matplotlib >= 2.0
* sklearn >= 0.20
```

with:

```markdown
## Dependency ##

ZFC requires Python 3.11 or newer and installs its runtime dependencies
(numpy, scipy, pandas, matplotlib, scikit-learn) automatically via `uv sync`
or `pip install`. Development dependencies (pytest, ruff) are installed with
`uv sync` from the `dev` dependency group.
```

- [ ] **Step 2: Replace the stale help block (lines 49-89)**

Replace the `usage: zfc [-h] [-i INPUT] ...` block up to the closing
```` ``` ```` with the current help (include `-i` as required):

```text
usage: zfc [-h] -i INPUT [-o OUTPREFIX]
           [--normalization {total,median,upper_quartile,median_ratio,none}]
           [--top-n-sgrna TOP_N_SGRNA] [--top-n-gene TOP_N_GENE]
           [--null-iteration NULL_ITERATION] [--seed SEED] [--plot] [--version]

Calculate fold change of screening data (zscore log2 fold change).

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Raw count table with header should include: <gene>,
                        <guide>, <barcode>, <ctrl>, <exp>.<ctrl> is the raw
                        counts of control group, <exp> is the raw counts of
                        treatment group. For screening without barcode, the
                        barcode column can be the same with guide.
  -o OUTPREFIX, --outprefix OUTPREFIX
                        Output file prefix, can be the file directory path
                        with file name prefix. The directory in the outprefix
                        should be built before analysis.
  --normalization {total,median,upper_quartile,median_ratio,none}
                        Normalization of raw count data, default is total.
                        Support method: total (Total sum normalization);
                        median (Median normalization); upper_quartile (Upper
                        quartile normalization (0.75)); median_ratio (Median
                        ratio normalization); none (No normalization).
  --top-n-sgrna TOP_N_SGRNA
                        Only consider top n barcodes for each sgRNA. Default
                        to use all the data.
  --top-n-gene TOP_N_GENE
                        Only consider top n barcodes for each gene. Default to
                        use all the data.
  --null-iteration NULL_ITERATION
                        The iteration to generate null distribution in
                        calculating the p value of genes. The larger the
                        better, but slower in calculation, default to be 100.
  --seed SEED            Random seed for the null-distribution permutations.
                        Default 0 (fixed) for reproducible output; use -1 for
                        nondeterministic output.
  --plot                 Output figures.
  --version              show program's version number and exit
```

- [ ] **Step 3: Verify the help text matches the binary**

Run: `uv run zfc --help`
Expected: the printed usage/options match the block above (note the leading
`-i INPUT` without brackets, reflecting the new `required=True`).

- [ ] **Step 4: Commit**

```bash
git add README.md
git commit -m "docs: sync README help and dependency sections"
```

---

## Definition of Done (run once at the end)

```bash
uvx ruff check
uv run pytest -q
```
Expected: ruff clean and all tests green (existing + new), including
`tests/test_parity.py::test_golden_output_parity`. Then confirm:

```bash
git status --short
git log --oneline -9
```
Expected: a clean tree and the nine commits above, in order.