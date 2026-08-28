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


def _lfc_std_model(bar_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, linear_model.LinearRegression]:
    """Fit the lfc-std linear model and add lfc_std (returns bar_df, train, reg)."""
    lfc_q1 = bar_df['lfc'].quantile(0.05)
    lfc_q2 = bar_df['lfc'].quantile(0.95)
    ctrl_q1 = bar_df['ctrl'].quantile(0.025)
    ctrl_q2 = bar_df['ctrl'].quantile(0.95)

    model_data = bar_df.loc[
        (bar_df['lfc'] >= lfc_q1) & (bar_df['lfc'] <= lfc_q2) &
        (bar_df['ctrl'] >= ctrl_q1) & (bar_df['ctrl'] <= ctrl_q2)
    ]

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


def _rank_aggregation(bar_df: pd.DataFrame, level: str, top_n: int | None) -> dict[str, pd.Series]:
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
    for rank_col, ascending in [('rank_down', True), ('rank_up', True)]:
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
                seed: int = 0) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, linear_model.LinearRegression]:
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