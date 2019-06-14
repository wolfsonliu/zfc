####################
# Statistics Functions
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################
import numpy as np
import functools


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


def df_smallcount(df, quantile=0.05, drop0=True):
    dfquantile = df.quantile(quantile, axis=0)
    smallidx = (df <= dfquantile).all(axis=1)
    smallcount = np.asarray(df[smallidx]).ravel()
    if drop0:
        smallcount = smallcount[smallcount > 0]
    return smallcount
