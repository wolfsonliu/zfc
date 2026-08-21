import numpy as np
import pandas as pd
import pytest

from zfc import zfoldchange


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
