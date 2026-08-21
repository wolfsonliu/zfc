import numpy as np
import pandas as pd

from zfc.statsfunc import df_normalization


def test_documented_methods_accepted():
    df = pd.DataFrame({
        'ctrl': [10.0, 20.0, 30.0],
        'exp': [5.0, 15.0, 25.0],
    })
    for method in ['total', 'median', 'upper_quartile', 'median_ratio', 'none']:
        result = df_normalization(df, method)
        assert result.shape == df.shape
        assert list(result.columns) == ['ctrl', 'exp']


def test_none_is_identity():
    df = pd.DataFrame({
        'ctrl': [10.0, 20.0, 30.0],
        'exp': [5.0, 15.0, 25.0],
    })
    assert df_normalization(df, 'none').equals(df)


def test_total_equalizes_column_sums():
    df = pd.DataFrame({
        'ctrl': [10.0, 20.0, 30.0],
        'exp': [5.0, 15.0, 25.0],
    })
    sums = df_normalization(df, 'total').sum(axis=0).to_numpy()
    assert np.isclose(sums[0], sums[1])