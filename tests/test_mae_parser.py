from pathlib import Path
import pytest
from openfep.mae_parser import parse_mae, CTBlock

FIXTURE = Path(__file__).parent / "fixtures" / "two_ct.mae"


def test_parse_returns_only_tagged_cts():
    blocks = parse_mae(FIXTURE)
    # Third CT has no s_fep_struc_tag, must be skipped
    assert len(blocks) == 2


def test_first_ct_is_receptor():
    blocks = parse_mae(FIXTURE)
    assert blocks[0].struc_tag == "receptor"


def test_second_ct_is_ligand():
    blocks = parse_mae(FIXTURE)
    assert blocks[1].struc_tag == "ligand"


def test_raw_block_contains_ct_header():
    blocks = parse_mae(FIXTURE)
    assert "f_m_ct" in blocks[0].raw_block


def test_ct_block_dataclass():
    b = CTBlock(struc_tag="receptor", raw_block="f_m_ct {}")
    assert b.struc_tag == "receptor"
    assert b.raw_block == "f_m_ct {}"


def test_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        parse_mae("/nonexistent/path.mae")
