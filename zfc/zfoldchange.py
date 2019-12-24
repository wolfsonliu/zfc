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
from .statsfunc import ecdf
from .statsfunc import p_adjust
from .statsfunc import df_robust_rank_aggregation
from .statsfunc import df_mean_rank_aggregation


def zfoldchange(data,
                top_n_sgrna=None,
                top_n_gene=None,
                iteration=100,
                normalization='total'):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in data.columns, 'data should have column: {}'.format(x)
    data.replace(np.nan, 0, inplace=True)
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

    # get valid data range
    lfc_q1 = bar_df['lfc'].quantile(0.05)
    lfc_q2 = bar_df['lfc'].quantile(0.95)
    ctrl_q1 = bar_df['ctrl'].quantile(0.025)
    ctrl_q2 = bar_df['ctrl'].quantile(0.95)

    # select valid data for the model
    model_data = bar_df.loc[
        (bar_df['lfc'] >= lfc_q1) & (bar_df['lfc'] <= lfc_q2) &
        (bar_df['ctrl'] >= ctrl_q1) & (bar_df['ctrl'] <= ctrl_q2)
    ]

    # cutting bins of the dat by ctrl data
    bins = pd.cut(model_data['ctrl'], 200).astype(str)

    train_data = pd.DataFrame(
        {
            'lfcstd': model_data.groupby(bins)['lfc'].std(),
            'ctrlmean': model_data.groupby(bins)['ctrl'].mean()
        }
    )

    train_data = train_data.loc[~train_data['lfcstd'].isnull()]

    # linear model of the lfc_std:ctrl_mean model
    reg = linear_model.LinearRegression()
    reg.fit(
        np.expand_dims(train_data['ctrlmean'].values, 1),
        train_data['lfcstd'].values
    )

    # calculate the lfc_std by ctrl_mean using parameters from the linear model
    bar_df = bar_df.assign(
        lfc_std=bar_df['ctrl'] * reg.coef_[0] + reg.intercept_
    )

    # adjust the lfc_stds that are less than or equal to zero
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

    # Calculate P value using normal distribution
    bar_df['p'] = stats.norm.cdf(bar_df['zlfc'])
    bar_df['p'] = bar_df['p'].map(
        lambda x: x if x <= 0.5 else 1 - x
    )

    # ------------------
    # Step 5: Calculate sgRNA mean zscore of fold change

    sg_df = pd.DataFrame(
        {
            'zlfc': bar_df.groupby(level=[0, 1])['zlfc'].mean(),
            'count': bar_df.groupby(level=[0, 1])['zlfc'].count()
        }
    )

    # generate null distribution of sgRNA zlfc
    sg_null_list = dict()
    sg_ecdf_list = dict()

    for i in range(iteration):
        sample_bar = bar_df[['zlfc']].copy()
        sample_bar.loc[:, 'zlfc'] = bar_df['zlfc'].sample(
            frac=1
        ).values
        sg_zlfc = sample_bar.groupby(
            level=[0, 1]
        ).mean()

        for c in sg_df['count'].unique():
            if c not in sg_null_list:
                sg_null_list[c] = list()
            sg_null_list[c].extend(
                sg_zlfc.loc[sg_df['count'] == c, 'zlfc'].values.tolist()
            )
        for c in sg_null_list:
            sg_ecdf_list[c] = ecdf(sg_null_list[c])

    # Calculate p value and FDR of sgRNA zlfc
    sg_df['tp'] = sg_df.apply(
        lambda a: sg_ecdf_list[a['count']](a['zlfc']),
        axis=1
    )
    sg_df['p'] = sg_df['tp'].map(
        lambda a: a if a <= 0.5 else 1 - a
    )
    sg_df['p_adj'] = p_adjust(sg_df['p'], 'BH')
    del sg_df['tp']
    del sg_null_list
    del sg_ecdf_list

    # -------------------
    # Step 6: Calculate Gene mean zscore of fold change

    g_df = pd.DataFrame(
        {
            'zlfc': sg_df.groupby(level=0)['zlfc'].mean(),
            'count': sg_df.groupby(level=0)['zlfc'].count(),
        }
    )

    # generate null distribution of gene zlfc
    g_null_list = dict()
    g_ecdf_list = dict()

    for i in range(iteration):
        sample_sg = sg_df[['zlfc']].copy()
        sample_sg.loc[:, 'zlfc'] = sg_df['zlfc'].sample(
            frac=1
        ).values
        g_zlfc = sample_sg.groupby(
            level=0
        ).mean()

        for c in g_df['count'].unique():
            if c not in g_null_list:
                g_null_list[c] = list()
            g_null_list[c].extend(
                g_zlfc.loc[g_df['count'] == c, 'zlfc'].values.tolist()
            )
        for c in g_null_list:
            g_ecdf_list[c] = ecdf(g_null_list[c])

    # Calculate p value and FDR of gene zlfc
    g_df['tp'] = g_df.apply(
        lambda a: g_ecdf_list[a['count']](a['zlfc']),
        axis=1
    )
    g_df['p'] = g_df['tp'].map(
        lambda a: a if a <= 0.5 else 1 - a
    )
    g_df['p_adj'] = p_adjust(g_df['p'], 'BH')
    del g_df['tp']
    del g_null_list
    del g_ecdf_list

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
