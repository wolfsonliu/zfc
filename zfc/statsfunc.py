####################
# Statistics Functions
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################
import pandas as pd
import numpy as np
import functools
import scipy.stats as stats


def ecdf(x):
    '''Calculate the empirical cumulative density function
       ====================
       Return a cumulative density function based
       on the input data.

    '''
    x = np.array(x)
    x.sort()
    n = len(x)
    y = np.concatenate([[0], np.linspace(1/n, 1, n)])

    def childfunc(sample, x, y, sorted=True):
        if not sorted:
            asort = np.argsort(x)
            x = np.take(x, asort, 0)
            y = np.take(y, asort, 0)
        idx = np.searchsorted(x, sample)
        return y[idx]

    return functools.partial(childfunc, x=x, y=y)


def bonferroni_correction(pvals):
    '''Bonferroni multitest p value correction'''

    pvals = np.asarray(pvals)
    n = len(pvals)
    pvals_corrected = pvals * float(n)
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected


def bh_correction(pvals):
    '''BH multitest p value correction'''

    pvals = np.asarray(pvals)
    n = len(pvals)
    pvals_decreasing_argsort = pvals.argsort()[::-1]
    pvals_reverse_argsort = pvals_decreasing_argsort.argsort()
    pvals_corrected = np.minimum.accumulate(
        pvals[pvals_decreasing_argsort] * n / np.arange(n, 0, -1)
    )[pvals_reverse_argsort]
    pvals_corrected[pvals_corrected > 1] = 1

    return pvals_corrected


def p_adjust(pvals, method='BH'):
    assert method in ['BH', 'Bonferroni'], 'Adjust method should be in: BH, Bonferroni'
    adjust_func = {
        'BH': bh_correction,
        'Bonferroni': bonferroni_correction
    }

    return adjust_func[method](pvals)


def df_median_ratio_normfactor(df):
    '''
    Normalize input data indicated by the label.
    Now only median ratio normalization is available.
    '''
    dfgm = np.exp(np.log(df + 1.0).sum(axis=1) / df.shape[1]) - 1.0
    dfgm[dfgm <= 0] = 1
    meanfactor = df.div(dfgm, axis=0)
    normfactor = 1 / meanfactor.median(axis=0)
    return normfactor


def df_total_count_normfactor(df):
    colsum = df.sum(axis=0)
    colsummean = colsum.sum() / df.shape[1]
    normfactor = colsummean / colsum
    return normfactor


def df_normalization(df, method):
    normfactor = df_total_count_normfactor(df)
    if method == 'none':
        normfactor = np.array([1]*len(normfactor))
    elif method == 'median':
        medianfactor = df_median_ratio_normfactor(df)
        if (medianfactor == 0).any():
            exit('Median factor is zero, using total count normalization')
        elif ((df == 0).sum(axis=0) / df.shape[0] > 0.45).any():
            exit('Too many zeros in counts, using total count normalization')
        else:
            normfactor = medianfactor
    elif method == 'total':
        pass
    result = df.mul(normfactor, axis=1)
    return result


def df_smallcount(df, quantile=0.1, drop0=True):
    dfquantile = df.quantile(quantile, axis=0)
    smallidx = (df <= dfquantile).all(axis=1)
    smallcount = np.asarray(df[smallidx]).ravel()
    if drop0:
        smallcount = smallcount[smallcount > 0]
    return smallcount


def mean_rank_aggregation(rank_matrix):
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


def rank_order_aggregation(rank_matrix):
    # Stuart, J.M., Segal, E., Koller, D., and Kim, S.K. (2003). A
    # Gene-Coexpression Network for Global Discovery of Conserved
    # Genetic Modules. Science 302, 249.

    # Aerts, S., Lambrechts, D., Maity, S., Van Loo, P., Coessens, B.,
    # De Smet, F., Tranchevent, L.-C., De Moor, B., Marynen, P.,
    # Hassan, B., et al. (2006). Gene prioritization through genomic
    # data fusion. Nat Biotech 24, 537–544.

    def Q(r):
        factorial = np.vectorize(np.math.factorial)
        r = np.array(r)
        r = r[~np.isnan(r)]
        r.sort()
        rr = -1 * r[::-1]
        k = len(rr)
        v = list()
        v.append(1)
        irange = np.arange(1, k + 1)
        for i in range(k):
            nowv = -1 * (
                np.array(v)[::-1] /
                factorial(irange[0:i + 1]) *
                (rr[i] ** irange[0:i + 1])
            ).sum()
            v.append(nowv)
        theQ = np.math.factorial(k) * v[-1]
        return theQ
    return rank_matrix.apply(lambda x: Q(x), axis=1)


def robust_rank_aggregation(rank_matrix):
    # Kolde, R., Laur, S., Adler, P., and Vilo, J. (2012). Robust rank
    # aggregation for gene list integration and
    # meta-analysis. Bioinformatics 28, 573–580.

    def beta_score(r):
        r = np.array(r)
        r.sort()
        rr = r[~np.isnan(r)]
        k = len(rr)
        return stats.beta.cdf(rr, np.arange(1, k + 1), np.arange(k, 0, -1))

    def rho_score(r):
        r = np.array(r)
        k = sum(~np.isnan(r))
        return (beta_score(r) * k).clip(0, 1).min()

    return rank_matrix.apply(lambda x: rho_score(x), axis=1)
