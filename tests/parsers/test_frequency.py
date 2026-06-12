import pytest
from turbomole_analyzer.parsers.frequency import FrequencyParser


@pytest.fixture
def parser():
    return FrequencyParser()


def test_parse_aoforce_style(tmp_path, parser):
    f = tmp_path / "aoforce.out"
    f.write_text("zero point vibrational energy :    0.541203 Hartree\n")
    assert parser.parse(f) == pytest.approx(0.541203)


def test_parse_aoforce_hyphen_style(tmp_path, parser):
    f = tmp_path / "aoforce.out"
    f.write_text("zero-point vibrational energy:    0.541203 Hartree\n")
    assert parser.parse(f) == pytest.approx(0.541203)


def test_parse_vibrational_zpe_style(tmp_path, parser):
    f = tmp_path / "aoforce.out"
    f.write_text("vibrational zero point energy :    0.123456 Hartree\n")
    assert parser.parse(f) == pytest.approx(0.123456)


def test_parse_job_last_style(tmp_path, parser):
    f = tmp_path / "job.last"
    f.write_text("Zero Point Energy :    0.0750000 Hartree\n")
    assert parser.parse(f) == pytest.approx(0.075)


def test_parse_case_insensitive(tmp_path, parser):
    f = tmp_path / "aoforce.out"
    f.write_text("ZERO POINT VIBRATIONAL ENERGY :    0.541203 Hartree\n")
    assert parser.parse(f) == pytest.approx(0.541203)


def test_parse_missing_file(tmp_path, parser):
    assert parser.parse(tmp_path / "aoforce.out") is None


def test_parse_no_zpe_data(tmp_path, parser):
    f = tmp_path / "aoforce.out"
    f.write_text("SCF energy:   -100.500000\nSome other output\n")
    assert parser.parse(f) is None
