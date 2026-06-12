import pytest
from turbomole_analyzer.parsers.nmr import NMRParser


@pytest.fixture
def parser():
    return NMRParser()


def test_parse_mpshift_format(tmp_path, parser):
    f = tmp_path / "mpshift.out"
    f.write_text("ATOM  C    1      ISOTROPIC:      120.500000       ANISOTROPIC:       50.000000\n")
    result = parser.parse(f)
    assert result is not None
    assert "C" in result
    assert result["C"]["1"] == pytest.approx(120.5)


def test_parse_job_last_format(tmp_path, parser):
    f = tmp_path / "job.last"
    f.write_text("Atom   1 C   isotropic shielding =    120.0 ppm\n")
    result = parser.parse(f)
    assert result is not None
    assert "C" in result
    assert result["C"]["1"] == pytest.approx(120.0)


def test_parse_negative_shielding(tmp_path, parser):
    f = tmp_path / "mpshift.out"
    f.write_text("ATOM  Pd   1      ISOTROPIC:      -277.672910       ANISOTROPIC:       338.359117\n")
    result = parser.parse(f)
    assert result["Pd"]["1"] == pytest.approx(-277.672910)


def test_parse_multiple_elements(tmp_path, parser):
    f = tmp_path / "mpshift.out"
    f.write_text(
        "ATOM  C    1      ISOTROPIC:      120.500000       ANISOTROPIC:       50.000000\n"
        "ATOM  H    2      ISOTROPIC:       31.450000       ANISOTROPIC:       10.000000\n"
    )
    result = parser.parse(f)
    assert "C" in result
    assert "H" in result
    assert result["H"]["2"] == pytest.approx(31.45)


def test_parse_multiple_atoms_same_element(tmp_path, parser):
    f = tmp_path / "mpshift.out"
    f.write_text(
        "ATOM  H    1      ISOTROPIC:       31.450000       ANISOTROPIC:       10.000000\n"
        "ATOM  H    2      ISOTROPIC:       30.120000       ANISOTROPIC:        9.800000\n"
    )
    result = parser.parse(f)
    assert len(result["H"]) == 2
    assert result["H"]["1"] == pytest.approx(31.45)
    assert result["H"]["2"] == pytest.approx(30.12)


def test_parse_missing_file(tmp_path, parser):
    assert parser.parse(tmp_path / "mpshift.out") is None


def test_parse_no_nmr_data(tmp_path, parser):
    f = tmp_path / "mpshift.out"
    f.write_text("Some output without NMR data\n")
    assert parser.parse(f) is None
