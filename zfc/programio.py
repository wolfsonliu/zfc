import pandas as pd


def read_raw_count(filepath):
    data = pd.read_csv(filepath, header=0, sep='\t')
