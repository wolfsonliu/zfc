import numpy as np
import pandas as pd
import pytest
from scipy.stats import beta

from zfc.statsfunc import df_robust_rank_aggregation
from zfc.statsfunc import p_adjust


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


def test_p_adjust_rejects_unknown_method():
    with pytest.raises(ValueError):
        p_adjust(np.array([0.01, 0.5]), 'bogus')


def test_p_adjust_returns_corrected_array():
    p = np.array([0.01, 0.05, 0.5])
    for method in ('BH', 'Bonferroni'):
        adj = p_adjust(p, method)
        assert adj.shape == p.shape
        assert ((adj >= 0) & (adj <= 1)).all()