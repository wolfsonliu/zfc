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
        'barcode_same_direction', 'lfc_std_modified', 'zlfc', 'p'
    ]].sort_values(
        ['gene', 'guide', 'barcode']
    ).to_csv(filepath, index=False, sep='\t')


def write_gresult(data, filepath):
    data[[
        'zlfc', 'count', 'p', 'FDR', 
        'RRA_Score_down', 'RRA_Score_up', 
        'MeanRank_Score_down', 'MeanRank_Score_up' 
    ]].sort_values(
        'zlfc', ascending=True
    ).to_csv(filepath, index=True, sep='\t')
