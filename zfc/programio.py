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
    return data
