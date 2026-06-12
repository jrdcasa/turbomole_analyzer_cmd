import pytest
from turbomole_analyzer.parsers.coord import CoordParser

BOHR_TO_ANG = 0.5291772109


@pytest.fixture
def parser():
    return CoordParser()


def test_parse_single_atom(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n 0.0 0.0 0.0 c\n$end")
    result = parser.parse(f)
    assert result is not None
    assert result.num_atoms == 1
    assert result.atoms[0].element == "C"


def test_parse_converts_bohr_to_angstrom(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n 1.0 2.0 3.0 h\n$end")
    result = parser.parse(f)
    assert result.atoms[0].x == pytest.approx(1.0 * BOHR_TO_ANG)
    assert result.atoms[0].y == pytest.approx(2.0 * BOHR_TO_ANG)
    assert result.atoms[0].z == pytest.approx(3.0 * BOHR_TO_ANG)


def test_parse_multiple_atoms(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n 0.0 0.0 0.0 c\n 1.889726 0.0 0.0 h\n$end")
    result = parser.parse(f)
    assert result.num_atoms == 2
    assert result.atoms[0].element == "C"
    assert result.atoms[1].element == "H"


def test_parse_element_cleaned_of_digits(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n 0.0 0.0 0.0 c1\n$end")
    result = parser.parse(f)
    assert result.atoms[0].element == "C"


def test_parse_skips_inline_comments(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n 0.0 0.0 0.0 c  # carbon\n$end")
    result = parser.parse(f)
    assert result.atoms[0].element == "C"


def test_parse_missing_file(tmp_path, parser):
    assert parser.parse(tmp_path / "coord") is None


def test_parse_empty_coord_block(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n$end")
    assert parser.parse(f) is None


def test_num_atoms_consistent_with_atoms_list(tmp_path, parser):
    f = tmp_path / "coord"
    f.write_text("$coord\n 0.0 0.0 0.0 c\n 1.0 0.0 0.0 h\n 2.0 0.0 0.0 h\n$end")
    result = parser.parse(f)
    assert result.num_atoms == len(result.atoms) == 3
