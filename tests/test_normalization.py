import numpy as np
import pandas as pd
import pytest

from zfc.statsfunc import df_normalization, df_quantile_normfactor


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


def test_upper_quartile_uses_quantile_factor():
    df = pd.DataFrame({
        'ctrl': [10.0, 20.0, 30.0, 25.0],
        'exp': [5.0, 15.0, 25.0, 10.0],
    })
    result = df_normalization(df, 'upper_quartile')
    expected = df.mul(df_quantile_normfactor(df), axis=1)
    assert result.equals(expected)
    assert not result.equals(df_normalization(df, 'total'))


def test_median_ratio_falls_back_to_total_on_many_zeros():
    df = pd.DataFrame({
        'ctrl': [0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
        'exp': [0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    })
    with pytest.warns(UserWarning, match='Too many zeros'):
        result = df_normalization(df, 'median_ratio')
    assert result.equals(df_normalization(df, 'total'))
