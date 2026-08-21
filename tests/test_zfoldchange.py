import numpy as np
import pandas as pd

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