import pandas as pd
import numpy as np
import scipy.stats as stats
from .statsfunc import df_normalization
from .statsfunc import df_smallcount
from .statsfunc import ecdf
from .statsfunc import bonferroni_correction
from .statsfunc import mean_rank_aggregation
from .statsfunc import rank_order_aggregation
from .statsfunc import robust_rank_aggregation


def zfoldchange(data,
                punish_rate=0.5,
                zero_sd_n=1.2,
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
    for a in [0.05, 0.1]:
        smallcount = df_smallcount(norm)
        if len(smallcount) > 0:
            break
    if len(smallcount) == 0:
        smallcount = np.array(
            [data['ctrl'][data['ctrl'] > 0].quantile(0.05)]
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

    refnorm = sgdf['ctrl']
    # calculate 5 bins from 10 bins cutting mean and std
    # to avoid the extreme value in positive screening
    bins10 = pd.cut(sgdf['lfc'], bins=10)

    lfcmean = sgdf['lfc'].loc[
        bins10.isin(bins10.value_counts(ascending=False).iloc[0:3].index)
    ].mean()
    lfcstd = sgdf['lfc'].loc[
        bins10.isin(bins10.value_counts(ascending=False).iloc[0:3].index)
    ].std()

    sidx = (
        sgdf['lfc'] <= (lfcmean + 2.5 * lfcstd)
    ) & (
        sgdf['lfc'] >= (lfcmean - 2.5 * lfcstd)
    )

    bins = pd.cut(refnorm, 100).astype(str)
    bins_count = bins.value_counts()

    bins_std = sgdf['lfc'][
        np.asarray(sidx)
    ].groupby(bins[np.asarray(sidx)]).std()
    bins_std[bins_count < (bins_count.sum() / 400)] = np.nan

    bins_mean = sgdf['lfc'][
        np.asarray(sidx)
    ].groupby(bins[np.asarray(sidx)]).mean()

    sgdf = sgdf.assign(
        lfc_mean=np.array(
            bins_mean[np.asarray(bins)]
        ),
        lfc_std=np.array(
            bins_std[np.asarray(bins)]
        )
    )

    na_std = sgdf['lfc'][sgdf['lfc_std'].isna()].std()
    sgdf.loc[sgdf['lfc_std'].isna(), 'lfc_std'] = na_std

    # ------------------
    # Step 4: Considering barcode direction

    sgdf['upper'] = (
        sgdf['lfc'] >= (sgdf['lfc_std'] * zero_sd_n)
    )

    sgdf['lower'] = (
        sgdf['lfc'] <= (sgdf['lfc_std'] * (-zero_sd_n))
    )

    sgdf['zero'] = (
        sgdf['lfc'] >= (sgdf['lfc_std'] * (-zero_sd_n))
    ) & (
        sgdf['lfc'] <= (sgdf['lfc_std'] * zero_sd_n)
    )

    sgdf['barcode_same_direction'] = np.asarray(
        sgdf.groupby('guide')[
            'upper', 'lower', 'zero'
        ].all().any(axis=1)[sgdf['guide']]
    )

    sgdf['max_lfc_std'] = np.asarray(
        sgdf.groupby('guide')[
            'lfc_std'
        ].max()[sgdf['guide']]
    )

    sgdf['lfc_std_modified'] = (
        1 - sgdf['barcode_same_direction'].astype(float)
    ) * sgdf['max_lfc_std'] * punish_rate + sgdf['lfc_std']

    del sgdf['max_lfc_std']
    del sgdf['upper']
    del sgdf['lower']
    del sgdf['zero']

    # ------------------
    # Step 5: Calculate zscore of fold change

    sgdf['zlfc'] = sgdf['lfc'] / sgdf['lfc_std_modified']

    # Calculate P value using normal distribution
    sgdf['p'] = stats.norm.cdf(sgdf['zlfc'])
    sgdf['p'] = sgdf['p'].map(
        lambda x: x if x <= 0.5 else 1 - x
    )

    # ------------------
    # Step 6: Calculate gene mean zscore of fold change

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

    # ------------------
    # Step 7: Rank aggregation
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
    gdf['MeanRank_Score_down'] = mean_rank_aggregation(sgdown)
    gdf['MeanRank_Score_up'] = mean_rank_aggregation(sgup)

    return (sgdf, gdf)
