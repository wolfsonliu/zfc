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
    assert (tmp_path / 'out_gene_volcano.png').exists()