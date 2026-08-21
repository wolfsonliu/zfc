import pandas as pd

from zfc import (
    read_raw_count,
    write_bar_result,
    write_sg_result,
    write_g_result,
)


SG_COLUMNS = [
    'zlfc', 'count', 'p', 'p_adj',
    'RRA_Score_down', 'RRA_Score_down_adj',
    'RRA_Score_up', 'RRA_Score_up_adj',
    'Mean_Rank_down', 'Mean_Rank_up',
]


def test_read_raw_count(tmp_path):
    path = tmp_path / 'counts.tsv'
    path.write_text(
        'gene\tguide\tbarcode\tctrl\texp\n'
        'A1BG\tG1\tB1\t100\t50\n'
        'A1BG\tG1\tB2\t100\t60\n'
    )
    df = read_raw_count(str(path))
    assert list(df.columns) == ['gene', 'guide', 'barcode', 'ctrl', 'exp']
    assert len(df) == 2
    assert list(df.index) == ['A1BG_G1_B1', 'A1BG_G1_B2']


def test_write_bar_result(tmp_path):
    bar = pd.DataFrame({
        'ctrl': [1.0], 'exp': [2.0], 'fc': [2.0], 'lfc': [1.0],
        'lfc_std': [0.5], 'zlfc': [2.0], 'p': [0.01],
    }, index=['A1BG_G1_B1'])
    path = tmp_path / 'out_barcode.txt'
    write_bar_result(bar, str(path))
    text = path.read_text()
    for col in ['ctrl', 'exp', 'fc', 'lfc', 'lfc_std', 'zlfc', 'p']:
        assert col in text


def test_write_sg_result(tmp_path):
    sg = pd.DataFrame({c: [1] for c in SG_COLUMNS}, index=['A1BG_G1'])
    path = tmp_path / 'out_sgrna.txt'
    write_sg_result(sg, str(path))
    text = path.read_text()
    for col in SG_COLUMNS:
        assert col in text


def test_write_g_result(tmp_path):
    g = pd.DataFrame({c: [1] for c in SG_COLUMNS}, index=['A1BG'])
    path = tmp_path / 'out_gene.txt'
    write_g_result(g, str(path))
    text = path.read_text()
    for col in SG_COLUMNS:
        assert col in text