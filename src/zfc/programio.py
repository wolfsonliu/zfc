####################
# ZFC
# Author: Wolfson Liu
# Email: wolfsonliu@live.com
####################

import pandas as pd


def read_raw_count(filepath,
                   gene_colname='gene',
                   guide_colname='guide',
                   barcode_colname='barcode',
                   ctrl_colname='ctrl',
                   exp_colname='exp'):
    data = pd.read_csv(filepath, header=0, sep='\t')
    data = data[
        [
            gene_colname, guide_colname, barcode_colname,
            ctrl_colname, exp_colname
        ]
    ]
    data.columns = ['gene', 'guide', 'barcode', 'ctrl', 'exp']
    data.index = data[['gene', 'guide', 'barcode']].apply(
        lambda a: '_'.join(a), axis=1
    )
    return data


def write_bar_result(data, filepath):
    data[[
        'ctrl', 'exp', 'fc', 'lfc', 'lfc_std', 'zlfc', 'p'
    ]].sort_index(
    ).to_csv(filepath, index=True, sep='\t')


def write_sg_result(data, filepath):
    data[[
        'zlfc', 'count',
        'p', 'p_adj',
        'RRA_Score_down', 'RRA_Score_down_adj',
        'RRA_Score_up', 'RRA_Score_up_adj',
        'Mean_Rank_down', 'Mean_Rank_up',
    ]].sort_values(
        'zlfc', ascending=True
    ).to_csv(filepath, index=True, sep='\t')


def write_g_result(data, filepath):
    data[[
        'zlfc', 'count',
        'p', 'p_adj',
        'RRA_Score_down', 'RRA_Score_down_adj',
        'RRA_Score_up', 'RRA_Score_up_adj',
        'Mean_Rank_down', 'Mean_Rank_up',
    ]].sort_values(
        'zlfc', ascending=True
    ).to_csv(filepath, index=True, sep='\t')
