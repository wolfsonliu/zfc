# ZFC Code Quality & Robustness Hardening Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use subagent-driven-development (recommended) or executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply six small, output-preserving code-quality/robustness changes to `statsfunc.py`, `zfoldchange.py`, and `plot.py` per `docs/specs/2026-08-28-zfc-code-quality-design.md`, with the existing golden parity test proving the three result tables stay bit-for-bit identical.

**Architecture:** Targeted edits to three modules plus two new tests. The only behavioral change is `p_adjust` raising `ValueError` instead of `AssertionError` on an unknown method; everything else is behavior-preserving (dead-code removal, type annotations, and two equivalent refactors guarded by the parity test `tests/test_parity.py`).

**Tech Stack:** Python 3.12+ (repo `>=3.11`), pandas, numpy, scipy, scikit-learn, ruff, pytest.

**Note on commands:** The project's canonical commands are `uv run pytest` and `uvx ruff check`, but in this sandbox `uv`/`uvx` fail on a read-only `~/.cache/uv`. Use the project venv directly: `.venv/bin/ruff check` and `.venv/bin/python -m pytest`. They are equivalent to the AGENTS.md commands.

---

## File Structure

| File | Change |
| --- | --- |
| `src/zfc/statsfunc.py` | Remove dead `ecdf()` + `import functools`; `p_adjust` `assert` → `ValueError`. |
| `src/zfc/zfoldchange.py` | Add return type annotations; group `_lfc_std_model` on the `pd.cut` Categorical instead of `.astype(str)`. |
| `src/zfc/plot.py` | Use `_log10p` in `zlfc_rra_scatter` / `zlfc_rrafdr_scatter`; use `np.linspace` in `lfcstd_ctrlmean_figure`. |
| `tests/test_statsfunc.py` | Add `p_adjust` tests (unknown method raises; known methods return corrected arrays). |
| `tests/test_zfoldchange.py` | Add `train_data` equivalence test proving the categorical grouping matches the old `.astype(str)` grouping. |

---

## Task 1: `p_adjust` raises `ValueError` on unknown method

**Files:**
- Modify: `src/zfc/statsfunc.py:68-75`
- Test: `tests/test_statsfunc.py`

- [ ] **Step 1: Write the failing test**

Append to `tests/test_statsfunc.py` and add the two imports at the top.

Add `import pytest` after `from scipy.stats import beta` and change the zfc import line to also pull in `p_adjust`:

```python
import numpy as np
import pandas as pd
import pytest
from scipy.stats import beta

from zfc.statsfunc import df_robust_rank_aggregation
from zfc.statsfunc import p_adjust
```

Append at the end of the file:

```python
def test_p_adjust_rejects_unknown_method():
    with pytest.raises(ValueError):
        p_adjust(np.array([0.01, 0.5]), 'bogus')


def test_p_adjust_returns_corrected_array():
    p = np.array([0.01, 0.05, 0.5])
    for method in ('BH', 'Bonferroni'):
        adj = p_adjust(p, method)
        assert adj.shape == p.shape
        assert ((adj >= 0) & (adj <= 1)).all()
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `.venv/bin/python -m pytest tests/test_statsfunc.py::test_p_adjust_rejects_unknown_method -v`
Expected: FAIL, with `AssertionError` raised (the current `assert` path), not `ValueError`.

- [ ] **Step 3: Implement `ValueError`**

In `src/zfc/statsfunc.py`, replace lines 68-69:

```python
def p_adjust(pvals, method: str = 'BH'):
    assert method in ['BH', 'Bonferroni'], 'Adjust method should be in: BH, Bonferroni'
```

with:

```python
def p_adjust(pvals, method: str = 'BH'):
    if method not in ('BH', 'Bonferroni'):
        raise ValueError('Adjust method should be in: BH, Bonferroni')
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `.venv/bin/python -m pytest tests/test_statsfunc.py -v`
Expected: all tests in the file PASS (the two new ones and the existing RRA test).

- [ ] **Step 5: Commit**

```bash
git add src/zfc/statsfunc.py tests/test_statsfunc.py
git commit -m "fix(statsfunc): p_adjust raises ValueError on unknown method"
```

---

## Task 2: Remove dead `ecdf()` and `import functools`

**Files:**
- Modify: `src/zfc/statsfunc.py:11,16-36`

No test exists for `ecdf` (verified: zero callers in `src/` and `tests/`); this is a dead-code removal guarded by the full suite.

- [ ] **Step 1: Remove the function and its import**

In `src/zfc/statsfunc.py`, delete the `import functools` line (line 11) and the entire `ecdf` function (lines 16-36) plus its trailing blank line, so the block that read:

```python
import pandas as pd
import numpy as np
import functools
import scipy.stats as stats
import warnings


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
```

now reads:

```python
import pandas as pd
import numpy as np
import scipy.stats as stats
import warnings


def bonferroni_correction(pvals: np.ndarray) -> np.ndarray:
```

- [ ] **Step 2: Verify nothing references `ecdf`**

Run: `grep -rn "ecdf" src/ tests/`
Expected: no matches (the inline `np.searchsorted` in `zfoldchange.py` uses variable names `ecdf_vals`, not the function — leave those).

- [ ] **Step 3: Run lint and the statsfunc tests**

Run: `.venv/bin/ruff check`
Expected: `All checks passed!`

Run: `.venv/bin/python -m pytest tests/test_statsfunc.py -v`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/zfc/statsfunc.py
git commit -m "refactor(statsfunc): remove dead ecdf function"
```

---

## Task 3: Complete type annotations in `zfoldchange.py`

**Files:**
- Modify: `src/zfc/zfoldchange.py:90,191,232-237`

No behavior change; verified by ruff + the full suite.

- [ ] **Step 1: Add return annotations**

Edit `src/zfc/zfoldchange.py`:

`_lfc_std_model` (line 90):

```python
def _lfc_std_model(bar_df: pd.DataFrame):
```

becomes:

```python
def _lfc_std_model(bar_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, linear_model.LinearRegression]:
```

`_rank_aggregation` (line 191):

```python
def _rank_aggregation(bar_df: pd.DataFrame, level: str, top_n: int | None) -> dict:
```

becomes:

```python
def _rank_aggregation(bar_df: pd.DataFrame, level: str, top_n: int | None) -> dict[str, pd.Series]:
```

`zfoldchange` (lines 232-237):

```python
def zfoldchange(data: pd.DataFrame,
                top_n_sgrna: int | None = None,
                top_n_gene: int | None = None,
                iteration: int = 100,
                normalization: str = 'total',
                seed: int = 0):
```

becomes:

```python
def zfoldchange(data: pd.DataFrame,
                top_n_sgrna: int | None = None,
                top_n_gene: int | None = None,
                iteration: int = 100,
                normalization: str = 'total',
                seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, linear_model.LinearRegression]:
```

`_empirical_p` (line 154) already has `-> pd.DataFrame`; leave its `groups` argument unannotated (no stable `SeriesGroupBy` type), as decided in the spec.

- [ ] **Step 2: Run lint and import smoke check**

Run: `.venv/bin/ruff check`
Expected: `All checks passed!`

Run: `.venv/bin/python -c "import zfc; print(zfc.__version__)"`
Expected: `0.2.2`

- [ ] **Step 3: Commit**

```bash
git add src/zfc/zfoldchange.py
git commit -m "refactor(zfoldchange): complete type annotations"
```

---

## Task 4: Group `_lfc_std_model` on the `pd.cut` Categorical

**Files:**
- Modify: `src/zfc/zfoldchange.py:102-111`
- Test: `tests/test_zfoldchange.py`

- [ ] **Step 1: Write the characterization test (passes on current code)**

In `tests/test_zfoldchange.py`, add the import:

```python
from zfc.zfoldchange import _lfc_std_model
```

and append:

```python
def _reference_lfc_std_train(bar_df):
    """Old `.astype(str)` grouping, kept as an independent reference."""
    lfc_q1 = bar_df['lfc'].quantile(0.05)
    lfc_q2 = bar_df['lfc'].quantile(0.95)
    ctrl_q1 = bar_df['ctrl'].quantile(0.025)
    ctrl_q2 = bar_df['ctrl'].quantile(0.95)
    model_data = bar_df.loc[
        (bar_df['lfc'] >= lfc_q1) & (bar_df['lfc'] <= lfc_q2) &
        (bar_df['ctrl'] >= ctrl_q1) & (bar_df['ctrl'] <= ctrl_q2)
    ]
    bins = pd.cut(model_data['ctrl'], 200).astype(str)
    train_data = pd.DataFrame({
        'lfcstd': model_data.groupby(bins)['lfc'].std(),
        'ctrlmean': model_data.groupby(bins)['ctrl'].mean(),
    })
    return train_data.loc[~train_data['lfcstd'].isnull()]


def test_lfc_std_model_grouping_unchanged():
    rng = np.random.default_rng(7)
    n = 8000
    bar = pd.DataFrame({
        'ctrl': rng.uniform(10.0, 4000.0, n),
        'lfc': rng.normal(0.0, 1.0, n),
    })

    _, train_new, _ = _lfc_std_model(bar)
    reference = _reference_lfc_std_train(bar)

    assert list(train_new.index) == list(reference.index)
    assert np.array_equal(train_new['ctrlmean'].to_numpy(), reference['ctrlmean'].to_numpy())
    assert np.array_equal(train_new['lfcstd'].to_numpy(), reference['lfcstd'].to_numpy())
```

- [ ] **Step 2: Run it to confirm it passes on the current code**

Run: `.venv/bin/python -m pytest tests/test_zfoldchange.py::test_lfc_std_model_grouping_unchanged -v`
Expected: PASS (current code equals the reference by construction).

- [ ] **Step 3: Implement the categorical grouping**

In `src/zfc/zfoldchange.py`, replace lines 102-111:

```python
    bins = pd.cut(model_data['ctrl'], 200).astype(str)

    train_data = pd.DataFrame(
        {
            'lfcstd': model_data.groupby(bins)['lfc'].std(),
            'ctrlmean': model_data.groupby(bins)['ctrl'].mean()
        }
    )

    train_data = train_data.loc[~train_data['lfcstd'].isnull()]
```

with:

```python
    bins = pd.cut(model_data['ctrl'], 200)

    train_data = pd.DataFrame(
        {
            'lfcstd': model_data.groupby(bins, observed=True)['lfc'].std(),
            'ctrlmean': model_data.groupby(bins, observed=True)['ctrl'].mean()
        }
    )
    train_data.index = train_data.index.astype(str)
    train_data = train_data.sort_index()

    train_data = train_data.loc[~train_data['lfcstd'].isnull()]
```

The `.sort_index()` is **required**: the old `groupby(astype(str))` sorted the
string bin labels lexicographically, while the categorical groupby sorts them
numerically. Without it the (order-sensitive) regression fit changes at the
last bit and the parity test fails.

- [ ] **Step 4: Run the characterization test and the golden parity test**

Run: `.venv/bin/python -m pytest tests/test_zfoldchange.py::test_lfc_std_model_grouping_unchanged tests/test_parity.py -v`
Expected: PASS (equivalence holds; the three result tables remain bit-for-bit identical).

- [ ] **Step 5: Commit**

```bash
git add src/zfc/zfoldchange.py tests/test_zfoldchange.py
git commit -m "perf(zfoldchange): group lfc-std model on categorical, not row-wise str"
```

---

## Task 5: Plot fixes (`_log10p` + `np.linspace`)

**Files:**
- Modify: `src/zfc/plot.py:234,261,271`

No new test (the plot smoke test `tests/test_plot.py` exercises these builders); guarded by ruff + `test_plot.py`.

- [ ] **Step 1: Use `_log10p` for the RRA scatters**

In `src/zfc/plot.py`, in `zlfc_rra_scatter` (line 261):

```python
    axes.scatter(data['zlfc'], np.log10(rra) * -1, color=BLUE, alpha=0.15, edgecolor='none')
```

becomes:

```python
    axes.scatter(data['zlfc'], _log10p(rra), color=BLUE, alpha=0.15, edgecolor='none')
```

In `zlfc_rrafdr_scatter` (line 271):

```python
    axes.scatter(data['zlfc'], np.log10(rra) * -1, color=BLUE, alpha=0.15, edgecolor='none')
```

becomes:

```python
    axes.scatter(data['zlfc'], _log10p(rra), color=BLUE, alpha=0.15, edgecolor='none')
```

- [ ] **Step 2: Use `np.linspace` for the regression line**

In `src/zfc/plot.py`, `lfcstd_ctrlmean_figure` (line 234):

```python
    x = np.arange(train_data['ctrlmean'].min(), train_data['ctrlmean'].max())
```

becomes:

```python
    x = np.linspace(train_data['ctrlmean'].min(), train_data['ctrlmean'].max(), 100)
```

- [ ] **Step 3: Run lint and the plot smoke test**

Run: `.venv/bin/ruff check`
Expected: `All checks passed!`

Run: `.venv/bin/python -m pytest tests/test_plot.py -v`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/zfc/plot.py
git commit -m "fix(plot): avoid log10(0)=inf in RRA scatters and degenerate std line"
```

---

## Task 6: Full verification

**Files:** none (verification only)

- [ ] **Step 1: Run the full lint + test suite**

Run: `.venv/bin/ruff check`
Expected: `All checks passed!`

Run: `.venv/bin/python -m pytest -q --no-header -p no:cacheprovider`
Expected: `63 passed` (or `63 passed` + the new tests, e.g. `66 passed`), no failures.

- [ ] **Step 2: Confirm parity explicitly**

Run: `.venv/bin/python -m pytest tests/test_parity.py -v`
Expected: PASS — barcode/sgRNA/gene tables are bit-for-bit identical to the `0.2.2` baseline.

- [ ] **Step 3: Final commit is not needed** (all changes are committed in Tasks 1-5); if any whitespace/trailing fix arose, commit it with `git commit -am "chore: final cleanup"`.