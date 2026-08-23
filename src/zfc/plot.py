####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .statsfunc import gini

# Blue theme palette (see docs/specs/2026-08-22-plot-theme-design.md).
BLUE = '#007BFF'
BLUE_DEEP = '#004999'
BLUE_MID = '#66ADFF'
BLUE_LIGHT = '#E8F1FE'
ACCENT = '#FF8C00'
GRID = '#D6DEE8'
TEXT_MUTED = '#5A6B7B'
NEUTRAL = '#9AA8B8'

# Default significance thresholds / limits for the result figures.
VOLCANO_FDR = 0.05
VOLCANO_X = 2.0
TOP_N_LABEL = 15
TOP_N_HITS = 20


def _apply_theme():
    plt.rcParams.update({
        'figure.dpi': 150,
        'savefig.dpi': 150,
        'font.size': 11,
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'axes.edgecolor': BLUE_DEEP,
        'axes.labelcolor': BLUE_DEEP,
        'axes.titlecolor': BLUE_DEEP,
        'xtick.color': TEXT_MUTED,
        'ytick.color': TEXT_MUTED,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'axes.axisbelow': True,
    })


def _style_ax(ax, grid=True):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if grid:
        ax.grid(axis='y', color=GRID, linewidth=0.5, alpha=0.6)


def _log10p(values):
    return -np.log10(np.maximum(values, np.finfo(float).tiny))


def _save(fig, outprefix, name):
    fig.tight_layout()
    fig.savefig('_'.join([outprefix, name + '.png']))
    fig.savefig('_'.join([outprefix, name + '.pdf']))
    plt.close(fig)


_apply_theme()


def lorenz_curve(data):
    countdata = dict()
    countdata['ctrl'] = data['ctrl'].sort_values(ascending=True).copy()
    countdata['exp'] = data['exp'].sort_values(ascending=True).copy()
    lorenz = dict()
    lorenz['ctrl'] = countdata['ctrl'].cumsum() / countdata['ctrl'].sum()
    lorenz['exp'] = countdata['exp'].cumsum() / countdata['exp'].sum()

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 6.4))

    curve_colors = {'ctrl': BLUE_DEEP, 'exp': BLUE_MID}
    for i, lab in enumerate(['ctrl', 'exp']):
        axes[i].plot(
            np.arange(lorenz[lab].size + 1) / lorenz[lab].size,
            np.insert(lorenz[lab].values, 0, 0),
            color=curve_colors[lab],
            linewidth=2
        )
        axes[i].plot([0, 1], [0, 1], color=ACCENT, linestyle='dotted', linewidth=1.5)
        axes[i].set_title(
            'Gini index: {0}'.format(round(gini(countdata[lab].values), 4))
        )
        axes[i].set_xlabel(['Control', 'Experiment'][i])
        axes[i].set_ylabel('Total count ratio')
        axes[i].set_aspect('equal')
        _style_ax(axes[i])

    return fig


def counts_boxplot(inputdata, barresult):
    fig, axes = plt.subplots(1, 1)
    box = axes.boxplot(
        [
            inputdata['ctrl'], inputdata['exp'],
            barresult['ctrl'], barresult['exp']
        ],
        patch_artist=True,
        sym='.'
    )
    for patch in box['boxes']:
        patch.set_facecolor(BLUE_LIGHT)
        patch.set_edgecolor(BLUE_DEEP)
        patch.set_linewidth(1.2)
    for part in ('whiskers', 'caps', 'medians'):
        for line in box[part]:
            line.set_color(BLUE_DEEP)
    for flier in box['fliers']:
        flier.set_markerfacecolor(BLUE)
        flier.set_markeredgecolor(BLUE)
        flier.set_markersize(3)
        flier.set_alpha(0.4)

    axes.set_title('Boxplot of raw counts and normalized counts')
    axes.set_xticklabels(['Raw Ctrl', 'Raw Exp', 'Norm Ctrl', 'Norm Exp'])
    _style_ax(axes)
    return fig


def normcount_scatter(data):
    fig, axes = plt.subplots(
        2, 2,
        gridspec_kw={'width_ratios': [3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(6.4, 6.4)
    )
    axes[1, 0].scatter(data['ctrl'], data['exp'], color=BLUE, alpha=0.15, edgecolor='none')
    xlim = axes[1, 0].get_xlim()
    ylim = axes[1, 0].get_ylim()
    lim = [min(xlim + ylim), max(xlim + ylim)]
    axes[1, 0].set_xlim(lim)
    axes[1, 0].set_ylim(lim)
    axes[1, 0].plot(lim, lim, linestyle='dotted', color=ACCENT, linewidth=1.5)
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel('Normalized counts of exp')

    axes[0, 0].hist(
        data['ctrl'], bins=100,
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6
    )
    axes[0, 0].set_ylabel('Counts')

    axes[1, 1].hist(
        data['exp'], bins=100,
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6,
        orientation='horizontal'
    )
    axes[1, 1].set_xlabel('Counts')

    axes[0, 1].set_visible(False)
    for ax in (axes[1, 0], axes[0, 0], axes[1, 1]):
        _style_ax(ax)
    return fig


def _lfc_scatter(data, value_col, ylabel):
    fig, axes = plt.subplots(
        2, 4,
        gridspec_kw={'width_ratios': [3, 1, 3, 1], 'height_ratios': [1, 3]},
        sharex='col',
        sharey='row',
        figsize=(12.8, 6.4)
    )
    axes[1, 0].scatter(data['ctrl'], data[value_col], color=BLUE, alpha=0.15, edgecolor='none')
    axes[1, 0].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 0].set_xlabel('Normalized counts of ctrl')
    axes[1, 0].set_ylabel(ylabel)
    axes[0, 0].hist(data['ctrl'], bins=100, color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6)
    axes[0, 0].set_ylabel('Counts')
    axes[1, 1].hist(
        data[value_col], bins=100,
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6,
        orientation='horizontal'
    )
    axes[1, 1].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 1].set_xlabel('Counts')
    axes[0, 1].set_visible(False)

    axes[1, 2].scatter(data['exp'], data[value_col], color=BLUE, alpha=0.15, edgecolor='none')
    axes[1, 2].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 2].set_xlabel('Normalized counts of exp')
    axes[1, 2].set_ylabel(ylabel)
    axes[0, 2].hist(data['exp'], bins=100, color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6)
    axes[0, 2].set_ylabel('Counts')
    axes[1, 3].hist(
        data[value_col], bins=100,
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6,
        orientation='horizontal'
    )
    axes[1, 3].axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes[1, 3].set_xlabel('Counts')
    axes[0, 3].set_visible(False)

    for ax in (axes[1, 0], axes[0, 0], axes[1, 1],
               axes[1, 2], axes[0, 2], axes[1, 3]):
        _style_ax(ax)
    return fig


def lfc_normcount_scatter(data):
    return _lfc_scatter(data, 'lfc', '$log_{2}$ fold change')


def zlfc_normcount_scatter(data):
    return _lfc_scatter(data, 'zlfc', 'Z score of $log_{2}$ fold change')


def lfcstd_lfc_scatter(data):
    fig, axes = plt.subplots(1, 1)
    axes.scatter(data['lfc_std'], data['lfc'], color=BLUE, alpha=0.15, edgecolor='none')
    axes.axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes.set_xlabel('Standard deviation of $log_{2}$ fold change')
    axes.set_ylabel('$log_{2}$ fold change')
    _style_ax(axes)
    return fig


def lfcstd_ctrlmean_figure(train_data, lm):
    fig, axes = plt.subplots(1, 1)
    axes.scatter(
        train_data['ctrlmean'], train_data['lfcstd'],
        color=BLUE, alpha=0.6, edgecolor='none'
    )
    x = np.arange(train_data['ctrlmean'].min(), train_data['ctrlmean'].max())
    y = x * lm.coef_[0] + lm.intercept_
    axes.plot(x, y, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes.set_title(
        'Slope: {0}; Intercept: {1}.'.format(
            round(lm.coef_[0], 4), round(lm.intercept_, 4)
        )
    )
    axes.set_xlabel('Control normalized count (bin mean)')
    axes.set_ylabel('Standard deviation of $log_{2}$ fold change (bin std)')
    _style_ax(axes)
    return fig


def zlfc_hist(data):
    fig, axes = plt.subplots(1, 1)
    axes.hist(data['zlfc'], bins=100, color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6)
    axes.axvline(0, linestyle='dotted', color=ACCENT, linewidth=1.5)
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('Counts')
    _style_ax(axes)
    return fig


def zlfc_rra_scatter(data):
    rra = data[['RRA_Score_down', 'RRA_Score_up']].min(axis=1)
    fig, axes = plt.subplots(1, 1)
    axes.scatter(data['zlfc'], np.log10(rra) * -1, color=BLUE, alpha=0.15, edgecolor='none')
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('$-log_{10}$(Robust Rank Aggregation)')
    _style_ax(axes)
    return fig


def zlfc_rrafdr_scatter(data):
    rra = data[['RRA_Score_down_adj', 'RRA_Score_up_adj']].min(axis=1)
    fig, axes = plt.subplots(1, 1)
    axes.scatter(data['zlfc'], np.log10(rra) * -1, color=BLUE, alpha=0.15, edgecolor='none')
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('$-log_{10}$(Robust Rank Aggregation FDR)')
    _style_ax(axes)
    return fig


def zlfc_volcano(data, top_n=None):
    zlfc = data['zlfc'].to_numpy()
    p_adj = data['p_adj'].to_numpy()
    y = _log10p(p_adj)
    sig_up = (zlfc > VOLCANO_X) & (p_adj < VOLCANO_FDR)
    sig_down = (zlfc < -VOLCANO_X) & (p_adj < VOLCANO_FDR)

    colors = np.full(len(data), NEUTRAL)
    colors[sig_up] = ACCENT
    colors[sig_down] = BLUE

    fig, axes = plt.subplots(1, 1, figsize=(6.4, 6.4))
    axes.scatter(zlfc, y, c=colors, s=12, alpha=0.35, edgecolor='none')
    axes.axhline(
        -np.log10(VOLCANO_FDR), color=BLUE_DEEP,
        linestyle='dashed', linewidth=1
    )
    axes.axvline(VOLCANO_X, color=BLUE_DEEP, linestyle='dashed', linewidth=1)
    axes.axvline(-VOLCANO_X, color=BLUE_DEEP, linestyle='dashed', linewidth=1)
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_ylabel('$-log_{10}$(adjusted p-value)')

    if top_n:
        sig_mask = sig_up | sig_down
        order = np.argsort(y)[::-1]
        labeled = 0
        for i in order:
            if not sig_mask[i]:
                continue
            name = data.index[i]
            if not isinstance(name, str):
                name = '_'.join(map(str, name))
            axes.annotate(
                str(name), (zlfc[i], y[i]),
                fontsize=8, color=BLUE_DEEP,
                xytext=(3, 3), textcoords='offset points'
            )
            labeled += 1
            if labeled >= top_n:
                break

    _style_ax(axes)
    return fig


def ma_plot(gresult, barresult):
    abun = (barresult['ctrl'] + barresult['exp']).groupby(level=0).mean() / 2
    x = np.log2(abun.reindex(gresult.index).to_numpy())
    y = gresult['zlfc'].to_numpy()

    fig, axes = plt.subplots(1, 1)
    axes.scatter(x, y, color=BLUE, s=12, alpha=0.35, edgecolor='none')
    axes.axhline(0, color=ACCENT, linestyle='dotted', linewidth=1.5)
    axes.set_xlabel('$log_{2}$(mean normalized count)')
    axes.set_ylabel('Z score of $log_{2}$ fold change')
    _style_ax(axes)
    return fig


def top_hits_rank(data, top_n=TOP_N_HITS):
    up = data.nlargest(top_n, 'zlfc')
    down = data.nsmallest(top_n, 'zlfc')
    sel = pd.concat([down, up]).drop_duplicates().sort_values('zlfc')

    fig, axes = plt.subplots(1, 1, figsize=(6.4, 0.28 * len(sel) + 2))
    colors = np.where(sel['zlfc'].to_numpy() < 0, BLUE, ACCENT)
    axes.barh(
        range(len(sel)), sel['zlfc'], color=colors,
        height=0.7, edgecolor='none'
    )
    axes.axvline(0, color=BLUE_DEEP, linewidth=1)

    ticks = [
        '{}{}'.format(str(g), '*' if p < VOLCANO_FDR else '')
        for g, p in zip(sel.index, sel['p_adj'])
    ]
    axes.set_yticks(range(len(sel)))
    axes.set_yticklabels(ticks, fontsize=8)
    axes.set_xlabel('Z score of $log_{2}$ fold change')
    axes.set_title('Top {} depleted / enriched genes'.format(top_n))
    _style_ax(axes)
    return fig


def p_hist(data):
    fig, axes = plt.subplots(1, 1)
    axes.hist(
        data['p'], bins=100, range=(0, 0.5),
        color=BLUE_LIGHT, edgecolor=BLUE, linewidth=0.6
    )
    expected = len(data) / 100
    axes.axhline(expected, color=ACCENT, linestyle='dashed', linewidth=1.2)
    axes.text(
        0.98, 0.95, 'uniform expectation',
        transform=axes.transAxes, ha='right', va='top',
        fontsize=8, color=TEXT_MUTED
    )
    axes.set_xlim(0, 0.5)
    axes.set_xlabel('p-value')
    axes.set_ylabel('Counts')
    axes.set_title('Distribution of permutation p-values')
    _style_ax(axes)
    return fig


def zfc_plot(outprefix, inputdata,
             barresult, sgresult, gresult,
             train_data, lm):
    _save(lorenz_curve(inputdata), outprefix, 'lorenz')
    _save(counts_boxplot(inputdata, barresult), outprefix, 'counts_boxplot')
    _save(normcount_scatter(barresult), outprefix, 'normcount_scatter')
    _save(lfc_normcount_scatter(barresult), outprefix, 'lfc_normcount_scatter')
    _save(zlfc_normcount_scatter(barresult), outprefix, 'zlfc_normcount_scatter')
    _save(lfcstd_lfc_scatter(barresult), outprefix, 'lfc_std_scatter')
    _save(zlfc_hist(sgresult), outprefix, 'sgrna_zlfc_hist')
    _save(zlfc_hist(gresult), outprefix, 'gene_zlfc_hist')
    _save(zlfc_rra_scatter(sgresult), outprefix, 'sgrna_zlfc_rra_scatter')
    _save(zlfc_rra_scatter(gresult), outprefix, 'gene_zlfc_rra_scatter')
    _save(zlfc_rrafdr_scatter(sgresult), outprefix, 'sgrna_zlfc_rra_fdr_scatter')
    _save(zlfc_rrafdr_scatter(gresult), outprefix, 'gene_zlfc_rra_fdr_scatter')
    _save(zlfc_volcano(gresult, top_n=TOP_N_LABEL), outprefix, 'gene_volcano')
    _save(zlfc_volcano(sgresult, top_n=None), outprefix, 'sgrna_volcano')
    _save(ma_plot(gresult, barresult), outprefix, 'gene_ma')
    _save(top_hits_rank(gresult), outprefix, 'gene_top_hits')
    _save(p_hist(sgresult), outprefix, 'sgrna_p_hist')
    _save(p_hist(gresult), outprefix, 'gene_p_hist')
    _save(lfcstd_ctrlmean_figure(train_data, lm), outprefix, 'lfcstd_ctrlmean')