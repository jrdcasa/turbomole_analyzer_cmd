import pytest
from turbomole_analyzer.parsers.gradient import GradientParser

BOHR_TO_ANG = 0.5291772109

SINGLE_CYCLE = """\
$grad
  cycle =   1    SCF energy =   -40.123456
   0.000000000000      0.000000000000      0.000000000000    c
   1.889726124565      0.000000000000      0.000000000000    h
  -1.23456789D-05   5.67890123D-06  -9.87654321D-07
   1.23456789D-05  -5.67890123D-06   9.87654321D-07
$end
"""

TWO_CYCLES = """\
$grad
  cycle =   1    SCF energy =   -40.123456
   0.000000000000      0.000000000000      0.000000000000    c
   1.889726124565      0.000000000000      0.000000000000    h
  -1.23456789D-05   5.67890123D-06  -9.87654321D-07
   1.23456789D-05  -5.67890123D-06   9.87654321D-07
  cycle =   2    SCF energy =   -40.125000
   0.001000000000      0.000000000000      0.000000000000    c
   1.888726124565      0.000000000000      0.000000000000    h
  -1.00000000D-06   0.00000000D+00   0.00000000D+00
   1.00000000D-06   0.00000000D+00   0.00000000D+00
$end
"""


@pytest.fixture
def parser():
    return GradientParser()


def test_parse_single_cycle_frame_count(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(SINGLE_CYCLE)
    frames = parser.parse_trajectory(f)
    assert len(frames) == 1


def test_parse_single_cycle_atom_count(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(SINGLE_CYCLE)
    frames = parser.parse_trajectory(f)
    assert len(frames[0]) == 2


def test_parse_element_labels(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(SINGLE_CYCLE)
    frames = parser.parse_trajectory(f)
    assert frames[0][0]["element"] == "C"
    assert frames[0][1]["element"] == "H"


def test_parse_converts_bohr_to_angstrom(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(SINGLE_CYCLE)
    frames = parser.parse_trajectory(f)
    assert frames[0][1]["x"] == pytest.approx(1.889726124565 * BOHR_TO_ANG)


def test_parse_fortran_d_notation(tmp_path, parser):
    content = (
        "$grad\n"
        "  cycle =   1\n"
        "   1.234567890123D+00   0.000000000000D+00   0.000000000000D+00    c\n"
        "  -1.00000000D-05   0.00000000D+00   0.00000000D+00\n"
        "$end\n"
    )
    f = tmp_path / "gradient"
    f.write_text(content)
    frames = parser.parse_trajectory(f)
    assert frames[0][0]["x"] == pytest.approx(1.234567890123 * BOHR_TO_ANG)


def test_parse_two_cycles(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(TWO_CYCLES)
    frames = parser.parse_trajectory(f)
    assert len(frames) == 2


def test_parse_gradient_lines_skipped(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(SINGLE_CYCLE)
    frames = parser.parse_trajectory(f)
    # Only coordinate lines (4 tokens) are captured; gradient lines (3 tokens) are skipped
    assert len(frames[0]) == 2


def test_parse_missing_file(tmp_path, parser):
    assert parser.parse_trajectory(tmp_path / "gradient") == []


def test_parse_empty_file(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text("")
    assert parser.parse_trajectory(f) == []


def test_parse_via_base_interface(tmp_path, parser):
    f = tmp_path / "gradient"
    f.write_text(SINGLE_CYCLE)
    assert parser.parse(f) == parser.parse_trajectory(f)
