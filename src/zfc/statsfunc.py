####################
# Statistics Functions
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

from __future__ import annotations

import pandas as pd
import numpy as np
import scipy.stats as stats
import warnings


def bonferroni_correction(pvals: np.ndarray) -> np.ndarray:
    '''Bonferroni multitest p value correction'''

    pvals = np.asarray(pvals)
    n = len(pvals)
    pvals_corrected = pvals * float(n)
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected


def bh_correction(pvals: np.ndarray) -> np.ndarray:
    '''BH multitest p value correction.

    Benjamini, Y., and Hochberg, Y. (1995). J. Roy. Stat. Soc. B 57, 289-300.
    '''

    pvals = np.asarray(pvals)
    n = len(pvals)
    pvals_decreasing_argsort = pvals.argsort()[::-1]
    pvals_reverse_argsort = pvals_decreasing_argsort.argsort()
    pvals_corrected = np.minimum.accumulate(
        pvals[pvals_decreasing_argsort] * n / np.arange(n, 0, -1)
    )[pvals_reverse_argsort]
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected


def p_adjust(pvals, method: str = 'BH'):
    if method not in ('BH', 'Bonferroni'):
        raise ValueError('Adjust method should be in: BH, Bonferroni')
    adjust_func = {
        'BH': bh_correction,
        'Bonferroni': bonferroni_correction
    }

    return adjust_func[method](pvals)


def df_median_ratio_normfactor(df: pd.DataFrame) -> pd.Series:
    '''
    Normalize input data indicated by the label.
    Now only median ratio normalization is available.

    Reference: Maza, E., Frasse, P., Senin, P., Bouzayen, M., and Zouine, M. (2013). Comparison of normalization methods for differential gene expression analysis in RNA-Seq experiments: A matter of relative size of studied transcriptomes. Commun Integr Biol 6, e25849–e25849.
    '''
    dfgm = np.exp(np.log(df + 1.0).sum(axis=1) / df.shape[1]) - 1.0
    dfgm[dfgm <= 0] = 1
    meanfactor = df.div(dfgm, axis=0)
    normfactor = 1 / meanfactor.median(axis=0)
    return normfactor


def df_total_count_normfactor(df: pd.DataFrame) -> pd.Series:
    colsum = df.sum(axis=0)
    normfactor = colsum.mean() / colsum
    return normfactor


def df_median_normfactor(df: pd.DataFrame) -> pd.Series:
    colmedian = df.median(axis=0)
    normfactor = colmedian.mean() / colmedian
    return normfactor


def df_quantile_normfactor(df: pd.DataFrame, quantile: float = 0.75) -> pd.Series:
    '''
    From edgeR .calcFactorQuantile
    '''
    normfactor = df_total_count_normfactor(df)
    y = df.mul(normfactor, axis=1) / df.sum(axis=0).mean()
    return y.quantile(quantile, axis=0)


def df_normalization(df: pd.DataFrame, method: str) -> pd.DataFrame:
    '''
    Normalization of raw count data.

    Normalization method:
        * "total": Total sum normalization
        * "median": Median normalization
        * "upper_quartile": Upper quartile normalization (0.75)
        * "median_ratio": Median ratio normalization
        * "none: No normalization
    '''
    assert method in [
        'total', 'median', 'upper_quartile', 'median_ratio', 'none'
    ], "Normalization method not support {0}".format(method)

    normfactor = df_total_count_normfactor(df)
    if method == 'none':
        normfactor = np.array([1]*len(normfactor))
    elif method == 'median_ratio':
        medianfactor = df_median_ratio_normfactor(df)
        if (medianfactor == 0).any():
            warnings.warn(
                'Median factor is zero; using total count normalization'
            )
        elif ((df == 0).sum(axis=0) / df.shape[0] > 0.45).any():
            warnings.warn(
                'Too many zeros in counts; using total count normalization'
            )
        else:
            normfactor = medianfactor
    elif method == 'upper_quartile':
        normfactor = df_quantile_normfactor(df)
    elif method == 'median':
        normfactor = df_median_normfactor(df)
    elif method == 'total':
        pass
    result = df.mul(normfactor, axis=1)
    return result


def df_smallcount(df: pd.DataFrame, quantile: float = 0.1,
                  drop0: bool = True) -> np.ndarray:
    dfquantile = df.quantile(quantile, axis=0)
    smallidx = (df <= dfquantile).all(axis=1)
    smallcount = np.asarray(df[smallidx]).ravel()
    if drop0:
        smallcount = smallcount[smallcount > 0]
    return smallcount


def df_mean_rank_aggregation(rank_matrix: pd.DataFrame) -> pd.Series:
    # Kolde, R., Laur, S., Adler, P., and Vilo, J. (2012). Robust rank
    # aggregation for gene list integration and
    # meta-analysis. Bioinformatics 28, 573–580.

    rank_mean = rank_matrix.mean(axis=1, skipna=True)
    rank_num = rank_matrix.notna().sum(axis=1)
    rank_table = pd.DataFrame({'rank': rank_mean, 'n': rank_num})
    # Rank is uniform distribution, with var = 1/12(max - min)
    # Central limit theorem, var_new = var / n
    rank_table['sd'] = np.sqrt(1/12/rank_table['n'])
    rank_score = pd.Series(
        stats.norm.cdf(
            rank_table['rank'], loc=0.5, scale=rank_table['sd']
        ).clip(0, 1)
    )
    rank_score.index = rank_matrix.index
    return rank_score


def df_robust_rank_aggregation(rank_matrix: pd.DataFrame) -> pd.Series:
    # Kolde, R., Laur, S., Adler, P., and Vilo, J. (2012). Robust rank
    # aggregation for gene list integration and
    # meta-analysis. Bioinformatics 28, 573–580.
    vals = np.asarray(rank_matrix, dtype=float)
    sorted_vals = np.sort(vals, axis=1)          # NaN sorts to the end
    k = (~np.isnan(sorted_vals)).sum(axis=1)     # non-NaN count per row
    n_cols = sorted_vals.shape[1]

    a = np.clip(np.arange(1, n_cols + 1), 1, None)      # 1..p per column
    b = np.clip(k[:, None] - np.arange(n_cols), 1, None)  # k-j per row/col
    beta_cdf = stats.beta.cdf(sorted_vals, a, b)

    masked = np.where(np.isnan(sorted_vals), np.nan, beta_cdf)
    scores = np.nanmin(np.clip(k[:, None] * masked, 0, 1), axis=1)
    return pd.Series(scores, index=rank_matrix.index)


def gini(arr: np.ndarray) -> float:
    '''
    Using R library ineq, Gini function
    '''
    sorted_arr = np.sort(arr)
    n = arr.size
    G = np.dot(sorted_arr, np.arange(n) + 1)
    G = 2 * G / sorted_arr.sum() - (n + 1.)
    return G / n
