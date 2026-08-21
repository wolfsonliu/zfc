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