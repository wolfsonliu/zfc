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
