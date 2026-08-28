import numpy as np
import pandas as pd
import pytest

from zfc import zfoldchange
from zfc.zfoldchange import _lfc_std_model


SG_COLUMNS = [
    'zlfc', 'count', 'p', 'p_adj',
    'RRA_Score_down', 'RRA_Score_down_adj',
    'RRA_Score_up', 'RRA_Score_up_adj',
    'Mean_Rank_down', 'Mean_Rank_up',
]


def _make_input(n_genes=50, n_guides=8, n_barcodes=4, seed=42):
    rng = np.random.default_rng(seed)
    rows = []
    for g in range(n_genes):
        gene = f'gene{g}'
        for s in range(n_guides):
            guide = f'guide{g}_{s}'
            for b in range(n_barcodes):
                barcode = f'bar{g}_{s}_{b}'
                ctrl = float(rng.integers(50, 4000))
                exp = ctrl * float(rng.uniform(0.3, 2.0))
                rows.append((gene, guide, barcode, ctrl, exp))
    return pd.DataFrame(rows, columns=['gene', 'guide', 'barcode', 'ctrl', 'exp'])


def test_zfoldchange_end_to_end():
    data = _make_input()
    bar, sg, g, train_data, _ = zfoldchange(data, iteration=5)
    assert len(bar) == len(data)
    for col in ['ctrl', 'exp', 'fc', 'lfc', 'lfc_std', 'zlfc', 'p']:
        assert col in bar.columns
    for col in SG_COLUMNS:
        assert col in sg.columns
        assert col in g.columns
    assert len(sg) == data.groupby(['gene', 'guide']).ngroups
    assert len(g) == data['gene'].nunique()
    assert len(train_data) > 100


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


def test_missing_column_raises_value_error():
    data = _make_input()
    data = data.drop(columns=['barcode'])
    with pytest.raises(ValueError, match='barcode'):
        zfoldchange(data, iteration=5, seed=0)


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
