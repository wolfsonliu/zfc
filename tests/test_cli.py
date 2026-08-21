import numpy as np
import pytest

from zfc.cli import main


def test_help_exits_zero():
    with pytest.raises(SystemExit) as exc_info:
        main(['--help'])
    assert exc_info.value.code == 0


def test_version_exits_zero():
    with pytest.raises(SystemExit) as exc_info:
        main(['--version'])
    assert exc_info.value.code == 0


def _write_counts(path):
    rng = np.random.default_rng(0)
    lines = ['gene\tguide\tbarcode\tctrl\texp']
    for g in range(50):
        for s in range(8):
            for b in range(4):
                ctrl = float(rng.integers(50, 4000))
                exp = ctrl * float(rng.uniform(0.3, 2.0))
                lines.append(
                    'gene{0}\tguide{0}_{1}\tbar{0}_{1}_{2}\t{3}\t{4}'.format(
                        g, s, b, ctrl, exp
                    )
                )
    path.write_text('\n'.join(lines) + '\n')


def test_seed_flag_produces_reproducible_output(tmp_path):
    counts = tmp_path / 'counts.tsv'
    _write_counts(counts)
    main(['-i', str(counts), '--seed', '3', '--outprefix', str(tmp_path / 'a')])
    first = (tmp_path / 'a_sgrna.txt').read_text()

    main(['-i', str(counts), '--seed', '3', '--outprefix', str(tmp_path / 'b')])
    second = (tmp_path / 'b_sgrna.txt').read_text()

    assert first == second


def test_missing_input_errors(tmp_path):
    with pytest.raises(SystemExit) as exc_info:
        main(['-i', str(tmp_path / 'nope.tsv'), '--outprefix', str(tmp_path / 'o')])
    assert exc_info.value.code == 2