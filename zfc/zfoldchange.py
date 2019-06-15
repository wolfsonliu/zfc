import pandas as pd
import numpy as np
import scipy.stats as stats
from .statsfunc import df_normalization
from .statsfunc import df_smallcount
from .statsfunc import ecdf
from .statsfunc import bonferroni_correction


def zfoldchange(data,
                punish_rate=0.5,
                iteration=100):
    # The data DF should contain: [gene, guide, barcode, ctrl, exp]

    for x in ['gene', 'guide', 'barcode', 'ctrl', 'exp']:
        assert x in data.columns, 'data should have column: {}'.format(x)

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
    sidx = (
        sgdf['lfc'] <= (
            sgdf['lfc'].mean() + 2.5 * sgdf['lfc'].std()
        )
    ) & (
        sgdf['lfc'] >= (
            sgdf['lfc'].mean() - 2.5 * sgdf['lfc'].std()
        )
    )

    bins = pd.cut(refnorm, 200).astype(str)

    bins_std = sgdf['lfc'][
        np.asarray(sidx)
    ].groupby(bins[np.asarray(sidx)]).std()
    bins_mean = sgdf['lfc'][
        np.asarray(sidx)
    ].groupby(bins[np.asarray(sidx)]).mean()
    bins_mean[bins_std.isna()] = 0
    bins_std = bins_std.replace(
        np.nan, sgdf['lfc'][np.asarray(sidx)].std()
    )
    sgdf = sgdf.assign(
        lfc_mean=np.array(
            bins_mean[np.asarray(bins)]
        ),
        lfc_std=np.array(
            bins_std[np.asarray(bins)]
        )
    )

    # ------------------
    # Step 4: Considering barcode direction

    sgdf['upper_or_zero'] = (
        sgdf['lfc'] >= (sgdf['lfc_std'] * (-1.2))
    )

    sgdf['lower_or_zero'] = (
        sgdf['lfc'] <= (sgdf['lfc_std'] * 1.2)
    )

    sgdf['barcode_same_direction'] = np.asarray(
        sgdf.groupby('guide')[
            'upper_or_zero', 'lower_or_zero'
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
    del sgdf['upper_or_zero']
    del sgdf['lower_or_zero']

    # ------------------
    # Step 5: Calculate zscore of fold change

    sgdf['zlfc'] = (
        sgdf['lfc'] - sgdf['lfc_mean']
    ) / sgdf['lfc_std_modified']

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

    return (sgdf, gdf)
