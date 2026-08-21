####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

import pandas as pd
import numpy as np
import scipy.stats as stats
from sklearn import linear_model
from .statsfunc import df_normalization
from .statsfunc import df_smallcount
from .statsfunc import p_adjust
from .statsfunc import df_robust_rank_aggregation
from .statsfunc import df_mean_rank_aggregation


def zfoldchange(data,
                top_n_sgrna=None,
                top_n_gene=None,
                iteration=100,
                normalization='total',
                seed=0):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in data.columns, 'data should have column: {}'.format(x)

    data = data.copy()

    # Validate identifiers and counts.
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

    # ------------------
    # Step 1: Normalization of raw counts
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
    bar_df = pd.concat(
        [
            data[['gene', 'guide', 'barcode']],
            norm
        ],
        axis=1, sort=False
    )
    del norm

    # ------------------
    # Step 2: Calculate fold change
    bar_df['fc'] = bar_df['exp'] / bar_df['ctrl']
    bar_df['lfc'] = np.log2(bar_df['fc'])

    # ------------------
    # Step 3: Calculate fold change std
    lfc_q1 = bar_df['lfc'].quantile(0.05)
    lfc_q2 = bar_df['lfc'].quantile(0.95)
    ctrl_q1 = bar_df['ctrl'].quantile(0.025)
    ctrl_q2 = bar_df['ctrl'].quantile(0.95)

    model_data = bar_df.loc[
        (bar_df['lfc'] >= lfc_q1) & (bar_df['lfc'] <= lfc_q2) &
        (bar_df['ctrl'] >= ctrl_q1) & (bar_df['ctrl'] <= ctrl_q2)
    ]

    bins = pd.cut(model_data['ctrl'], 200).astype(str)

    train_data = pd.DataFrame(
        {
            'lfcstd': model_data.groupby(bins)['lfc'].std(),
            'ctrlmean': model_data.groupby(bins)['ctrl'].mean()
        }
    )

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

    # ------------------
    # Step 4: Calculate raw zlfc using lfc and lfc_std
    bar_df.loc[:, 'zlfc'] = bar_df['lfc'] / bar_df['lfc_std']

    bar_df['rank_down'] = bar_df['zlfc'].rank(
        method='average', ascending=True
    ) / len(bar_df['zlfc'])

    bar_df['rank_up'] = bar_df['zlfc'].rank(
        method='average', ascending=False
    ) / len(bar_df['zlfc'])

    bar_df.set_index(['gene', 'guide', 'barcode'], inplace=True)

    # Calculate P value using normal distribution (vectorized one-sided p).
    bar_p = stats.norm.cdf(bar_df['zlfc'])
    bar_df['p'] = np.minimum(bar_p, 1 - bar_p)

    rng = (
        np.random.default_rng()
        if (seed is None or seed < 0)
        else np.random.default_rng(seed)
    )

    # ------------------
    # Step 5: Calculate sgRNA mean zscore of fold change
    sg_groups = bar_df.groupby(level=[0, 1])['zlfc']
    sg_df = pd.DataFrame(
        {
            'zlfc': sg_groups.mean(),
            'count': sg_groups.count()
        }
    )
    sg_codes = sg_groups.ngroup().to_numpy()
    sg_counts = np.bincount(sg_codes)
    sg_uniq = np.unique(sg_counts)
    sg_nulls = {int(c): [] for c in sg_uniq}
    sg_vals = bar_df['zlfc'].to_numpy()

    for _ in range(iteration):
        perm = rng.permutation(sg_vals)
        means = np.bincount(
            sg_codes, weights=perm, minlength=len(sg_df)
        ) / sg_counts
        for c in sg_uniq:
            sg_nulls[int(c)].append(means[sg_counts == c])
    for c in sg_uniq:
        sg_nulls[int(c)] = np.sort(np.concatenate(sg_nulls[int(c)]))

    sg_z = sg_df['zlfc'].to_numpy()
    sg_p = np.empty(len(sg_df), dtype=float)
    for c in sg_uniq:
        mask = sg_counts == c
        null = sg_nulls[int(c)]
        ecdf_vals = np.searchsorted(null, sg_z[mask], side='left') / null.size
        sg_p[mask] = np.minimum(ecdf_vals, 1 - ecdf_vals)
    sg_df['p'] = sg_p
    sg_df['p_adj'] = p_adjust(sg_df['p'], 'BH')
    del sg_nulls

    # -------------------
    # Step 6: Calculate Gene mean zscore of fold change
    g_groups = sg_df.groupby(level=0)['zlfc']
    g_df = pd.DataFrame(
        {
            'zlfc': g_groups.mean(),
            'count': g_groups.count()
        }
    )
    g_codes = g_groups.ngroup().to_numpy()
    g_counts = np.bincount(g_codes)
    g_uniq = np.unique(g_counts)
    g_nulls = {int(c): [] for c in g_uniq}
    g_vals = sg_df['zlfc'].to_numpy()

    for _ in range(iteration):
        perm = rng.permutation(g_vals)
        means = np.bincount(
            g_codes, weights=perm, minlength=len(g_df)
        ) / g_counts
        for c in g_uniq:
            g_nulls[int(c)].append(means[g_counts == c])
    for c in g_uniq:
        g_nulls[int(c)] = np.sort(np.concatenate(g_nulls[int(c)]))

    g_z = g_df['zlfc'].to_numpy()
    g_p = np.empty(len(g_df), dtype=float)
    for c in g_uniq:
        mask = g_counts == c
        null = g_nulls[int(c)]
        ecdf_vals = np.searchsorted(null, g_z[mask], side='left') / null.size
        g_p[mask] = np.minimum(ecdf_vals, 1 - ecdf_vals)
    g_df['p'] = g_p
    g_df['p_adj'] = p_adjust(g_df['p'], 'BH')
    del g_nulls

    # ------------------
    # Step 7: Rank aggregation

    # sgRNA Rank
    if top_n_sgrna is None:
        top_n_sgrna = int(
            bar_df.groupby(level=[0, 1])['zlfc'].count().median()
        )
    sg_b_down = bar_df[['rank_down']].copy()
    sg_b_down.loc[:, 'groupid'] = sg_b_down.groupby(
        level=[0, 1]
    )['rank_down'].rank(
        method='first', ascending=True
    ).astype(int)
    sg_b_down.reset_index(drop=False, inplace=True)
    del sg_b_down['barcode']
    sg_b_down.set_index(['gene', 'guide', 'groupid'], inplace=True)
    sg_b_down = sg_b_down.unstack(level=2)
    sg_b_down.columns = sg_b_down.columns.levels[1]
    sg_b_down = sg_b_down[list(range(1, top_n_sgrna + 1))]

    sg_b_up = bar_df[['rank_up']].copy()
    sg_b_up.loc[:, 'groupid'] = sg_b_up.groupby(
        level=[0, 1]
    )['rank_up'].rank(
        method='first', ascending=True
    ).astype(int)
    sg_b_up.reset_index(drop=False, inplace=True)
    del sg_b_up['barcode']
    sg_b_up.set_index(['gene', 'guide', 'groupid'], inplace=True)
    sg_b_up = sg_b_up.unstack(level=2)
    sg_b_up.columns = sg_b_up.columns.levels[1]
    sg_b_up = sg_b_up[list(range(1, top_n_sgrna + 1))]

    sg_df.loc[:, 'RRA_Score_down'] = df_robust_rank_aggregation(sg_b_down)
    sg_df.loc[:, 'RRA_Score_down_adj'] = p_adjust(
        sg_df['RRA_Score_down'], 'BH'
    )
    sg_df.loc[:, 'RRA_Score_up'] = df_robust_rank_aggregation(sg_b_up)
    sg_df.loc[:, 'RRA_Score_up_adj'] = p_adjust(
        sg_df['RRA_Score_up'], 'BH'
    )

    sg_df.loc[:, 'Mean_Rank_down'] = df_mean_rank_aggregation(sg_b_down)
    sg_df.loc[:, 'Mean_Rank_up'] = df_mean_rank_aggregation(sg_b_up)

    # gene Rank
    if top_n_gene is None:
        top_n_gene = int(
            bar_df.groupby(level=0)['zlfc'].count().median()
        )
    g_b_down = bar_df[['rank_down']].copy()
    g_b_down.loc[:, 'groupid'] = g_b_down.groupby(
        level=0
    )['rank_down'].rank(
        method='first', ascending=True
    ).astype(int)
    g_b_down.reset_index(drop=False, inplace=True)
    del g_b_down['barcode']
    del g_b_down['guide']
    g_b_down.set_index(['gene', 'groupid'], inplace=True)
    g_b_down = g_b_down.unstack(level=1)
    g_b_down.columns = g_b_down.columns.levels[1]
    g_b_down = g_b_down[list(range(1, top_n_gene + 1))]

    g_b_up = bar_df[['rank_up']].copy()
    g_b_up.loc[:, 'groupid'] = g_b_up.groupby(
        level=0
    )['rank_up'].rank(
        method='first', ascending=True
    ).astype(int)
    g_b_up.reset_index(drop=False, inplace=True)
    del g_b_up['barcode']
    del g_b_up['guide']
    g_b_up.set_index(['gene', 'groupid'], inplace=True)
    g_b_up = g_b_up.unstack(level=1)
    g_b_up.columns = g_b_up.columns.levels[1]
    g_b_up = g_b_up[list(range(1, top_n_gene + 1))]

    g_df.loc[:, 'RRA_Score_down'] = df_robust_rank_aggregation(g_b_down)
    g_df.loc[:, 'RRA_Score_down_adj'] = p_adjust(
        g_df['RRA_Score_down'], 'BH'
    )
    g_df.loc[:, 'RRA_Score_up'] = df_robust_rank_aggregation(g_b_up)
    g_df.loc[:, 'RRA_Score_up_adj'] = p_adjust(
        g_df['RRA_Score_up'], 'BH'
    )

    g_df.loc[:, 'Mean_Rank_down'] = df_mean_rank_aggregation(g_b_down)
    g_df.loc[:, 'Mean_Rank_up'] = df_mean_rank_aggregation(g_b_up)

    return (bar_df, sg_df, g_df, train_data, reg)