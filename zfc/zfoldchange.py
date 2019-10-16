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
from .statsfunc import bonferroni_correction
from .statsfunc import robust_rank_aggregation


def zfoldchange(data,
                leverage_threshold=0.7,
                topn=None,
                iteration=100):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in data.columns, 'data should have column: {}'.format(x)
    data.replace(np.nan, 0, inplace=True)
    # ------------------
    # Step 1: Normalization of raw counts
    norm = df_normalization(
        data[['ctrl', 'exp']],
        'total'
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
    sgdf = pd.concat(
        [
            data[['gene', 'guide', 'barcode']],
            norm
        ],
        axis=1, sort=False
    )
    del norm
    # ------------------
    # Step 2: Calculate fold change
    sgdf['fc'] = sgdf['exp'] / sgdf['ctrl']
    sgdf['lfc'] = np.log2(sgdf['fc'])

    # ------------------
    # Step 3: Calculate fold change std

    # get valid data range
    lfc_q1 = sgdf['lfc'].quantile(0.05)
    lfc_q2 = sgdf['lfc'].quantile(0.95)
    ctrl_q1 = sgdf['ctrl'].quantile(0.025)
    ctrl_q2 = sgdf['ctrl'].quantile(0.95)

    # select valid data for the model
    model_data = sgdf.loc[
        (sgdf['lfc'] >= lfc_q1) & (sgdf['lfc'] <= lfc_q2) &
        (sgdf['ctrl'] >= ctrl_q1) & (sgdf['ctrl'] <= ctrl_q2)
    ]

    # cutting bins of the dat by ctrl data
    bins = pd.cut(model_data['ctrl'], 200).astype(str)

    bins_lfcstd = model_data.groupby(bins)['lfc'].std()
    bins_ctrlmean = model_data.groupby(bins)['ctrl'].mean()

    null_idx = bins_lfcstd.isnull().copy()
    bins_ctrlmean = bins_ctrlmean[~null_idx]
    bins_lfcstd = bins_lfcstd[~null_idx]

    # linear model of the lfc_std:ctrl_mean model
    reg = linear_model.LinearRegression()
    reg.fit(np.expand_dims(bins_ctrlmean.values, 1), bins_lfcstd.values)

    # calculate the lfc_std by ctrl_mean using parameters from the linear model
    sgdf = sgdf.assign(
        lfc_std=sgdf['ctrl'] * reg.coef_[0] + reg.intercept_
    )

    # adjust the lfc_stds that are less than or equal to zero
    sgdf.loc[sgdf['lfc_std'] <= 0, 'lfc_std'] = sgdf.loc[
        sgdf['lfc_std'] > 0,
        'lfc_std'
    ].median()

    # ------------------
    # Step 4: Calculate raw zlfc using lfc and lfc_std

    sgdf['raw_zlfc'] = sgdf['lfc'] / sgdf['lfc_std']

    # ------------------
    # Step 5: Modify sgRNA-iBAR ZLFC with large leverage

    sgdf['sgrna_zlfc_mean'] = np.asarray(
        sgdf.groupby('guide')['raw_zlfc'].mean()[sgdf['guide']]
    )
    sgdf['sgrna_zlfc_var'] = np.asarray(
        sgdf.groupby('guide')['raw_zlfc'].var()[sgdf['guide']]
    )
    sgdf['sgrna_zlfc_count'] = np.asarray(
        sgdf.groupby('guide')['raw_zlfc'].count()[sgdf['guide']]
    )
    # calculate barcode leverage
    sgdf['barcode_zlfc_leverage'] = (
        sgdf['raw_zlfc'] - sgdf['sgrna_zlfc_mean']
    ) ** 2 / (
        sgdf['sgrna_zlfc_var'] * (sgdf['sgrna_zlfc_count'] - 1)
    ) + (1 / sgdf['sgrna_zlfc_count'])

    sgdf['barcode_zlfc_leverage_ratio'] = sgdf['sgrna_zlfc_count'] * (
        2 - sgdf['barcode_zlfc_leverage']
    ) / (2 * sgdf['sgrna_zlfc_count'] - 2)

    sgdf = sgdf.loc[
        sgdf['barcode_zlfc_leverage'] <= leverage_threshold
    ]

    # calculate zlfc from raw_zlfc and leverage ratio
    sgdf['zlfc'] = sgdf['raw_zlfc'] * sgdf['barcode_zlfc_leverage_ratio']

    # replace the infinit number with null
    sgdf = sgdf.replace([np.inf, -np.inf], np.nan)

    # assign raw_zlfc to zlfc if the zlfc is null
    sgdf.loc[sgdf['zlfc'].isnull(), 'zlfc'] = sgdf.loc[
        sgdf['zlfc'].isnull(), 'raw_zlfc'
    ]

    # ------------------
    # Step 6: Calculate zscore of fold change p value in normal distribution

    # Calculate P value using normal distribution
    sgdf['p'] = stats.norm.cdf(sgdf['zlfc'])
    sgdf['p'] = sgdf['p'].map(
        lambda x: x if x <= 0.5 else 1 - x
    )

    # ------------------
    # Step 7: Calculate gene mean zscore of fold change

    gdf = pd.DataFrame(
        {
            'zlfc': sgdf.groupby('gene')['zlfc'].mean(),
            'count': sgdf.groupby('gene')['zlfc'].count(),
        }
    )

    # generate null distribution of gene zlfc
    null_data = dict()
    ecdf_list = dict()

    for i in range(iteration):
        zlfc = sgdf['zlfc'].sample(
            frac=1
        ).reset_index(
            drop=True
        ).groupby(
            sgdf['gene'].values
        ).mean()

        for c in gdf['count'].unique():
            if c not in null_data:
                null_data[c] = list()
            null_data[c].extend(
                list(zlfc[gdf['count'] == c])
            )
        for c in null_data:
            ecdf_list[c] = ecdf(null_data[c])

    # Calculate p value and FDR of gene zlfc
    gdf['tp'] = gdf.apply(
        lambda a: ecdf_list[a['count']](a['zlfc']),
        axis=1
    )
    gdf['p'] = gdf['tp'].map(
        lambda a: a if a <= 0.5 else 1 - a
    )
    gdf['FDR'] = bonferroni_correction(gdf['p'])
    del gdf['tp']

    gdf = gdf.assign(
        zlfc_up=sgdf.loc[sgdf['zlfc'] > 0].groupby('gene')['zlfc'].mean(),
        count_up=sgdf.loc[sgdf['zlfc'] > 0].groupby('gene')['zlfc'].count(),
        zlfc_down=sgdf.loc[sgdf['zlfc'] < 0].groupby('gene')['zlfc'].mean(),
        count_down=sgdf.loc[sgdf['zlfc'] < 0].groupby('gene')['zlfc'].count(),
    )

    # ------------------
    # Step 8: Rank aggregation
    if topn is None:
        topn = int(sgdf.groupby('gene')['barcode'].count().median())

    sgdf['rank_down'] = sgdf['zlfc'].rank(
        method='average', ascending=True
    ) / len(sgdf['zlfc'])
    sgdown = sgdf[['gene', 'rank_down']]
    sgdown = sgdown.assign(
        groupid=sgdown.groupby('gene')['rank_down'].rank(
            method='first', ascending=True
        ).astype(int)
    ).set_index(['gene', 'groupid'])
    sgdown = sgdown.unstack()
    sgdown.columns = sgdown.columns.levels[1]
    sgdown = sgdown[list(range(1, topn + 1))]

    sgdf['rank_up'] = sgdf['zlfc'].rank(
        method='average', ascending=False
    ) / len(sgdf['zlfc'])
    sgup = sgdf[['gene', 'rank_up']]
    sgup = sgup.assign(
        groupid=sgup.groupby('gene')['rank_up'].rank(
            method='first', ascending=False
        ).astype(int)
    ).set_index(['gene', 'groupid'])
    sgup = sgup.unstack()
    sgup.columns = sgup.columns.levels[1]
    sgup = sgup[list(range(1, topn + 1))]

    gdf['RRA_Score_down'] = robust_rank_aggregation(sgdown)
    gdf['RRA_Score_up'] = robust_rank_aggregation(sgup)

    return (sgdf, gdf)
