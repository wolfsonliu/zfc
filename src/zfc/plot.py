####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

import matplotlib 
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from .statsfunc import gini


def lorenz_curve(data):
    countdata = dict()
    countdata['ctrl'] = data['ctrl'].sort_values(ascending=True).copy()
    countdata['exp'] = data['exp'].sort_values(ascending=True).copy()
    lorenz = dict()
    lorenz['ctrl'] = countdata['ctrl'].cumsum() / countdata['ctrl'].sum()
    lorenz['exp'] = countdata['exp'].cumsum() / countdata['exp'].sum()

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))

    for i, lab in enumerate(['ctrl', 'exp']):
        axes[i].plot(
            np.arange(lorenz[lab].size + 1) / lorenz[lab].size,
            np.insert(lorenz[lab].values, 0, 0),
            color='k'
        )
        axes[i].plot([0, 1], [0, 1], color='red', linestyle='dotted')
        axes[i].set_title(
            'Gini index: {0}'.format(round(gini(countdata[lab].values), 4))
        )
        axes[i].set_xlabel(['Control', 'Experiment'][i])
        axes[i].set_ylabel('Total count ratio')
        axes[i].set_aspect('equal')

    return fig


def counts_boxplot(inputdata, sgresult):
    fig, axes = plt.subplots(1, 1)
    axes.boxplot(
        [
            inputdata['ctrl'], inputdata['exp'],
            sgresult['ctrl'], sgresult['exp']
        ],
        sym='.k'
    )
    axes.set_title('Boxplot of raw counts and normalized counts')
    axes.set_xticklabels(['Raw Ctrl', 'Raw Exp', 'Norm Ctrl', 'Norm Exp'])
    return fig


def normcount_scatter(data):
    fig, axes = plt.subplots(
        2, 2,
        gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(6.4, 6.4)
    )
    axes[1, 0].scatter(
        data['ctrl'], data['exp'],
        c=['black'] * len(data['ctrl']), alpha=0.1,
        edgecolor=['none'] * len(data['ctrl'])
    )
    xlim = axes[1, 0].get_xlim()
    ylim = axes[1, 0].get_ylim()
    lim = [min(xlim + ylim), max(xlim + ylim)]
    axes[1, 0].set_xlim(lim)
    axes[1, 0].set_ylim(lim)
    axes[1, 0].plot(lim, lim, linestyle='dotted', color='red')
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel('Normalized counts of exp')

    axes[0, 0].hist(
        data['ctrl'], bins=100,
        color='black'
    )
    axes[0, 0].set_ylabel('Counts')

    axes[1, 1].hist(
        data['exp'], bins=100,
        color='black', orientation='horizontal'
    )
    axes[1, 1].set_xlabel('Counts')

    axes[0, 1].set_visible(False)
    return fig


def lfc_normcount_scatter(data):
    fig, axes = plt.subplots(
        2, 4,
        gridspec_kw={'width_ratios': [3, 1, 3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(12.8, 6.4)
    )
    axes[1, 0].scatter(
        data['ctrl'], data['lfc'],
        c=['black'] * len(data['ctrl']), alpha=0.1,
        edgecolor=['none'] * len(data['ctrl'])
    )
    axes[1, 0].axhline(0, color='red', linestyle='dotted')
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel('$log_{2}$ fold change')
    axes[0, 0].hist(data['ctrl'], bins=100, color='black')
    axes[0, 0].set_ylabel('Counts')
    axes[1, 1].hist(
        data['lfc'], bins=100,
        color='black', orientation='horizontal'
    )
    axes[1, 1].axhline(0, color='red', linestyle='dotted')
    axes[1, 1].set_xlabel('Counts')
    axes[0, 1].set_visible(False)

    axes[1, 2].scatter(
        data['exp'], data['lfc'],
        c=['black'] * len(data['exp']), alpha=0.1,
        edgecolor=['none'] * len(data['exp'])
    )
    axes[1, 2].axhline(0, color='red', linestyle='dotted')
    axes[1, 2].set_xlabel('Normalized counts of exp')
    axes[1, 2].set_ylabel('$log_{2}$ fold change')
    axes[0, 2].hist(data['exp'], bins=100, color='black')
    axes[0, 2].set_ylabel('Counts')
    axes[1, 3].hist(
        data['lfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 3].axhline(0, color='red', linestyle='dotted')
    axes[1, 3].set_xlabel('Counts')
    axes[0, 3].set_visible(False)

    plt.tight_layout()
    return fig


def zlfc_normcount_scatter(data):
    fig, axes = plt.subplots(
        2, 4,
        gridspec_kw={'width_ratios': [3, 1, 3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(12.8, 6.4)
    )
    axes[1, 0].scatter(
        data['ctrl'], data['zlfc'],
        c=['black'] * len(data['ctrl']), alpha=0.1,
        edgecolor=['none'] * len(data['ctrl'])
    )
    axes[1, 0].axhline(0, color='red', linestyle='dotted')
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel('Z score of $log_{2}$ fold change')
    axes[0, 0].hist(data['ctrl'], bins=100, color='black')
    axes[0, 0].set_ylabel('Counts')
    axes[1, 1].hist(
        data['zlfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 1].axhline(0, color='red', linestyle='dotted')
    axes[1, 1].set_xlabel('Counts')
    axes[0, 1].set_visible(False)

    axes[1, 2].scatter(
        data['exp'], data['zlfc'],
        c=['black'] * len(data['exp']), alpha=0.1,
        edgecolor=['none'] * len(data['exp'])
    )
    axes[1, 2].axhline(0, color='red', linestyle='dotted')
    axes[1, 2].set_xlabel('Normalized counts of exp')
    axes[1, 2].set_ylabel('Z score of $log_{2}$ fold change')
    axes[0, 2].hist(data['exp'], bins=100, color='black')
    axes[0, 2].set_ylabel('Counts')
    axes[1, 3].hist(
        data['zlfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 3].axhline(0, color='red', linestyle='dotted')
    axes[1, 3].set_xlabel('Counts')
    axes[0, 3].set_visible(False)

    plt.tight_layout()
    return fig


def lfcstd_lfc_scatter(data):
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        data['lfc_std'], data['lfc'],
        c=['black'] * len(data['lfc_std']), alpha=0.1,
        edgecolor=['none'] * len(data['lfc_std'])
    )
    axes.axhline(0, color='red', linestyle='dotted')
    axes.set_xlabel('Standard deviation of $log_{2}$ fold change')
    axes.set_ylabel('$log_{2}$ fold change')
    return fig


def lfcstd_ctrlmean_figure(train_data, lm):
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        train_data['ctrlmean'], train_data['lfcstd'], color='k'
    )
    x = np.arange(train_data['ctrlmean'].min(), train_data['ctrlmean'].max())
    y = x * lm.coef_[0] + lm.intercept_
    axes.plot(x, y, color='red', linestyle='dotted')
    axes.set_title(
        'Slope: {0}; Intercept: {1}.'.format(
            round(lm.coef_[0], 4), round(lm.intercept_, 4)
        )
    )
    axes.set_xlabel('Control normalized count (bin mean)')
    axes.set_ylabel('Standard deviation of $log_{2}$ fold change (bin std)')
    return fig


def zlfc_hist(data):
    fig, axes = plt.subplots(1, 1)
    axes.hist(data['zlfc'], bins=100, color='black')
    axes.axvline(0, linestyle='dotted', color='red')
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('Counts')
    return fig


def zlfc_rra_scatter(data):
    data.loc[:, 'RRA'] = data[['RRA_Score_down', 'RRA_Score_up']].min(axis=1)
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        data['zlfc'], np.log10(data['RRA']) * -1,
        c=['black'] * len(data['zlfc']), alpha=0.1,
        edgecolor=['none'] * len(data['zlfc'])
    )
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('Robust Rank Aggregation (-log10)')
    return fig


def zlfc_rrafdr_scatter(data):
    data.loc[:, 'RRA'] = data[
        ['RRA_Score_down_adj', 'RRA_Score_up_adj']
    ].min(axis=1)
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        data['zlfc'], np.log10(data['RRA']) * -1,
        c=['black'] * len(data['zlfc']), alpha=0.1,
        edgecolor=['none'] * len(data['zlfc'])
    )
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('Robust Rank Aggregation (-log10 FDR)')
    return fig


def zfc_plot(outprefix, inputdata,
             barresult, sgresult, gresult,
             train_data, lm):
    # Lorenz Curve
    fig = lorenz_curve(inputdata)
    fig.savefig('_'.join([outprefix, 'lorenz.png']))
    fig.savefig('_'.join([outprefix, 'lorenz.pdf']))

    # Boxplot: Counts and normalized counts
    fig = counts_boxplot(inputdata, barresult)
    fig.savefig('_'.join([outprefix, 'counts_boxplot.png']))
    fig.savefig('_'.join([outprefix, 'counts_boxplot.pdf']))

    # Scatterplot: Normalized counts of control and experiment
    fig = normcount_scatter(barresult)
    fig.savefig('_'.join([outprefix, 'normcount_scatter.png']))
    fig.savefig('_'.join([outprefix, 'normcount_scatter.pdf']))

    # Scatterplot: Barcode LFC with control and experiment normalized counts
    fig = lfc_normcount_scatter(barresult)
    fig.savefig('_'.join([outprefix, 'lfc_normcount_scatter.png']))
    fig.savefig('_'.join([outprefix, 'lfc_normcount_scatter.pdf']))

    # Scatterplot: Barcode ZLFC with control and experiment normalized counts
    fig = zlfc_normcount_scatter(barresult)
    fig.savefig('_'.join([outprefix, 'zlfc_normcount_scatter.png']))
    fig.savefig('_'.join([outprefix, 'zlfc_normcount_scatter.pdf']))

    # Scatterplot: Barcode LFC std and LFC
    fig = lfcstd_lfc_scatter(barresult)
    fig.savefig('_'.join([outprefix, 'lfc_std_scatter.png']))
    fig.savefig('_'.join([outprefix, 'lfc_std_scatter.pdf']))

    # Histogram: sgRNA ZLFC
    fig = zlfc_hist(sgresult)
    fig.savefig('_'.join([outprefix, 'sgrna_zlfc_hist.png']))
    fig.savefig('_'.join([outprefix, 'sgrna_zlfc_hist.pdf']))

    # Histogram: Gene ZLFC
    fig = zlfc_hist(gresult)
    fig.savefig('_'.join([outprefix, 'gene_zlfc_hist.png']))
    fig.savefig('_'.join([outprefix, 'gene_zlfc_hist.pdf']))

    # Scatterplot: sgRNA ZLFC RRA
    fig = zlfc_rra_scatter(sgresult)
    fig.savefig('_'.join([outprefix, 'sgrna_zlfc_rra_scatter.png']))
    fig.savefig('_'.join([outprefix, 'sgrna_zlfc_rra_scatter.pdf']))

    # Scatterplot: Gene ZLFC RRA
    fig = zlfc_rra_scatter(gresult)
    fig.savefig('_'.join([outprefix, 'gene_zlfc_rra_scatter.png']))
    fig.savefig('_'.join([outprefix, 'gene_zlfc_rra_scatter.pdf']))

    # Scatterplot: sgRNA ZLFC RRA FDR
    fig = zlfc_rrafdr_scatter(sgresult)
    fig.savefig('_'.join([outprefix, 'sgrna_zlfc_rra_fdr_scatter.png']))
    fig.savefig('_'.join([outprefix, 'sgrna_zlfc_rra_fdr_scatter.pdf']))

    # Scatterplot: Gene ZLFC RRA
    fig = zlfc_rrafdr_scatter(gresult)
    fig.savefig('_'.join([outprefix, 'gene_zlfc_rra_fdr_scatter.png']))
    fig.savefig('_'.join([outprefix, 'gene_zlfc_rra_fdr_scatter.pdf']))

    # Model figure of lfcstd and ctrlmean
    fig = lfcstd_ctrlmean_figure(train_data, lm)
    fig.savefig('_'.join([outprefix, 'lfcstd_ctrlmean.png']))
    fig.savefig('_'.join([outprefix, 'lfcstd_ctrlmean.pdf']))
