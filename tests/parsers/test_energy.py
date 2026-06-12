import pytest
from turbomole_analyzer.parsers.energy import EnergyParser


@pytest.fixture
def parser():
    return EnergyParser()


def test_parse_single_cycle(tmp_path, parser):
    f = tmp_path / "energy"
    f.write_text("$energy\n 1  -100.500000\n$end")
    assert parser.parse(f) == pytest.approx(-100.500000)


def test_parse_returns_last_cycle(tmp_path, parser):
    f = tmp_path / "energy"
    f.write_text("$energy\n 1  -100.500000\n 2  -100.501000\n$end")
    assert parser.parse(f) == pytest.approx(-100.501000)


def test_parse_negative_energy(tmp_path, parser):
    f = tmp_path / "energy"
    f.write_text("$energy\n 1  -277.672910\n$end")
    assert parser.parse(f) == pytest.approx(-277.672910)


def test_parse_zero_energy(tmp_path, parser):
    f = tmp_path / "energy"
    f.write_text("$energy\n 1  0.000000\n$end")
    assert parser.parse(f) == pytest.approx(0.0)


def test_parse_missing_file(tmp_path, parser):
    assert parser.parse(tmp_path / "energy") is None


def test_parse_empty_block(tmp_path, parser):
    f = tmp_path / "energy"
    f.write_text("$energy\n$end")
    assert parser.parse(f) is None
