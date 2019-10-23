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


def write_sgresult(data, filepath):
    data[[
        'gene', 'guide', 'barcode', 'ctrl', 'exp', 'fc', 'lfc', 'lfc_std',
        'barcode_zlfc_leverage', 'zlfc', 'p'
    ]].sort_values(
        ['gene', 'guide', 'barcode']
    ).to_csv(filepath, index=False, sep='\t')


def write_gresult(data, filepath):
    data[[
        'zlfc', 'zlfc_down', 'zlfc_up',
        'count', 'count_up', 'count_down',
        'p', 'p_adj',
        'RRA_Score_down', 'RRA_Score_down_adj',
        'RRA_Score_up', 'RRA_Score_up_adj',
        'Mean_Rank_down',
        'Mean_Rank_up',
    ]].sort_values(
        'zlfc', ascending=True
    ).to_csv(filepath, index=True, sep='\t')
