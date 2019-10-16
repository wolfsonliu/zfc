import matplotlib.pyplot as plt


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


def normcount_scatter(sgresult):
    fig, axes = plt.subplots(
        2, 2,
        gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(6.4, 6.4)
    )
    axes[1, 0].scatter(
        sgresult['ctrl'], sgresult['exp'],
        color='black', alpha=0.1, edgecolor=''
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
        sgresult['ctrl'], bins=100, color='black'
    )
    axes[0, 0].set_ylabel('Counts')

    axes[1, 1].hist(
        sgresult['exp'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 1].set_xlabel('Counts')

    axes[0, 1].set_visible(False)
    return fig


def lfc_normcount_scatter(sgresult):
    fig, axes = plt.subplots(
        2, 4,
        gridspec_kw={'width_ratios': [3, 1, 3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(12.8, 6.4)
    )
    axes[1, 0].scatter(
        sgresult['ctrl'], sgresult['lfc'],
        color='black', alpha=0.1, edgecolor='')
    axes[1, 0].axhline(0, color='red', linestyle='dotted')
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel('$log_{2}$ fold change')
    axes[0, 0].hist(sgresult['ctrl'], bins=100, color='black')
    axes[0, 0].set_ylabel('Counts')
    axes[1, 1].hist(
        sgresult['lfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 1].axhline(0, color='red', linestyle='dotted')
    axes[1, 1].set_xlabel('Counts')
    axes[0, 1].set_visible(False)

    axes[1, 2].scatter(
        sgresult['exp'], sgresult['lfc'],
        color='black', alpha=0.1, edgecolor=''
    )
    axes[1, 2].axhline(0, color='red', linestyle='dotted')
    axes[1, 2].set_xlabel('Normalized counts of exp')
    axes[1, 2].set_ylabel('$log_{2}$ fold change')
    axes[0, 2].hist(sgresult['exp'], bins=100, color='black')
    axes[0, 2].set_ylabel('Counts')
    axes[1, 3].hist(
        sgresult['lfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 3].axhline(0, color='red', linestyle='dotted')
    axes[1, 3].set_xlabel('Counts')
    axes[0, 3].set_visible(False)

    plt.tight_layout()
    return fig


def zlfc_normcount_scatter(sgresult):
    fig, axes = plt.subplots(
        2, 4,
        gridspec_kw={'width_ratios': [3, 1, 3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(12.8, 6.4)
    )
    axes[1, 0].scatter(
        sgresult['ctrl'], sgresult['zlfc'],
        color='black', alpha=0.1, edgecolor='')
    axes[1, 0].axhline(0, color='red', linestyle='dotted')
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel('Zscore of $log_{2}$ fold change')
    axes[0, 0].hist(sgresult['ctrl'], bins=100, color='black')
    axes[0, 0].set_ylabel('Counts')
    axes[1, 1].hist(
        sgresult['zlfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 1].axhline(0, color='red', linestyle='dotted')
    axes[1, 1].set_xlabel('Counts')
    axes[0, 1].set_visible(False)

    axes[1, 2].scatter(
        sgresult['exp'], sgresult['zlfc'],
        color='black', alpha=0.1, edgecolor=''
    )
    axes[1, 2].axhline(0, color='red', linestyle='dotted')
    axes[1, 2].set_xlabel('Normalized counts of exp')
    axes[1, 2].set_ylabel('Zscore of $log_{2}$ fold change')
    axes[0, 2].hist(sgresult['exp'], bins=100, color='black')
    axes[0, 2].set_ylabel('Counts')
    axes[1, 3].hist(
        sgresult['zlfc'], bins=100, color='black', orientation='horizontal'
    )
    axes[1, 3].axhline(0, color='red', linestyle='dotted')
    axes[1, 3].set_xlabel('Counts')
    axes[0, 3].set_visible(False)

    plt.tight_layout()
    return fig


def lfc_std_scatter(sgresult):
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        sgresult['lfc_std'], sgresult['lfc'],
        color='black', alpha=0.1, edgecolor=''
    )
    axes.axhline(0, color='red', linestyle='dotted')
    axes.set_xlabel('Standard deviation of $log_{2}$ fold change')
    axes.set_ylabel('$log_{2}$ fold change')
    return fig


def gene_zlfc_hist(gresult):
    fig, axes = plt.subplots(1, 1)
    axes.hist(gresult['zlfc'], bins=100, color='black')
    axes.axvline(0, linestyle='dotted', color='red')
    axes.set_xlabel('Zscore of $log_{2}$ fold change')
    axes.set_ylabel('Counts')
    return fig


def zfc_plot(outprefix, inputdata, sgresult, gresult):
    fig = counts_boxplot(inputdata, sgresult)
    fig.savefig('_'.join([outprefix, 'counts_boxplot.png']))
    fig.savefig('_'.join([outprefix, 'counts_boxplot.pdf']))

    fig = normcount_scatter(sgresult)
    fig.savefig('_'.join([outprefix, 'normcount_scatter.png']))
    fig.savefig('_'.join([outprefix, 'normcount_scatter.pdf']))

    fig = lfc_normcount_scatter(sgresult)
    fig.savefig('_'.join([outprefix, 'lfc_normcount_scatter.png']))
    fig.savefig('_'.join([outprefix, 'lfc_normcount_scatter.pdf']))

    fig = zlfc_normcount_scatter(sgresult)
    fig.savefig('_'.join([outprefix, 'zlfc_normcount_scatter.png']))
    fig.savefig('_'.join([outprefix, 'zlfc_normcount_scatter.pdf']))

    fig = lfc_std_scatter(sgresult)
    fig.savefig('_'.join([outprefix, 'lfc_std_scatter.png']))
    fig.savefig('_'.join([outprefix, 'lfc_std_scatter.pdf']))

    fig = gene_zlfc_hist(gresult)
    fig.savefig('_'.join([outprefix, 'gene_zlfc_hist.png']))
    fig.savefig('_'.join([outprefix, 'gene_zlfc_hist.pdf']))
