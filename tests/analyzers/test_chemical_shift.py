import pytest
from pathlib import Path

from turbomole_analyzer.analyzers.chemical_shift import ChemicalShiftCalculator, parse_element_values
from turbomole_analyzer.analyzers.workflow import WorkflowAnalyzer
from turbomole_analyzer.models.results import JobResults, NMRData


# ---------------------------------------------------------------------------
# parse_element_values
# ---------------------------------------------------------------------------

def test_parse_single_element():
    result = parse_element_values("C=188.1")
    assert result == {"C": pytest.approx(188.1)}


def test_parse_multiple_elements():
    result = parse_element_values("C=188.1,H=31.7")
    assert result["C"] == pytest.approx(188.1)
    assert result["H"] == pytest.approx(31.7)


def test_parse_normalises_element_case():
    result = parse_element_values("c=188.1,h=31.7")
    assert "C" in result
    assert "H" in result


def test_parse_ignores_whitespace():
    result = parse_element_values(" C = 188.1 , H = 31.7 ")
    assert result["C"] == pytest.approx(188.1)
    assert result["H"] == pytest.approx(31.7)


def test_parse_negative_value():
    result = parse_element_values("Pd=-277.67")
    assert result["Pd"] == pytest.approx(-277.67)


def test_parse_integer_value():
    result = parse_element_values("C=188")
    assert result["C"] == pytest.approx(188.0)


def test_parse_invalid_format_raises():
    with pytest.raises(ValueError):
        parse_element_values("C188.1")


def test_parse_invalid_value_raises():
    with pytest.raises(ValueError):
        parse_element_values("C=abc")


# ---------------------------------------------------------------------------
# ChemicalShiftCalculator — fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def nmr_data():
    return NMRData(
        method="nmr",
        chemical_shifts={
            "C": {"1": 120.5, "2": 135.0},
            "H": {"3": 28.3, "4": 29.1},
            "N": {"5": 100.0},
        },
    )


@pytest.fixture
def calc_c_only():
    return ChemicalShiftCalculator(sigma_ref={"C": 188.1})


# ---------------------------------------------------------------------------
# ChemicalShiftCalculator — basic formula
# ---------------------------------------------------------------------------

def test_calculate_default_delta_ref_zero(nmr_data, calc_c_only):
    result = calc_c_only.calculate(nmr_data)
    assert result["C"]["1"] == pytest.approx(188.1 - 120.5)
    assert result["C"]["2"] == pytest.approx(188.1 - 135.0)


def test_calculate_with_custom_delta_ref(nmr_data):
    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1}, delta_ref={"C": 5.0})
    result = calc.calculate(nmr_data)
    assert result["C"]["1"] == pytest.approx(5.0 + 188.1 - 120.5)


def test_calculate_multiple_elements(nmr_data):
    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1, "H": 31.7})
    result = calc.calculate(nmr_data)
    assert "C" in result
    assert "H" in result
    assert result["H"]["3"] == pytest.approx(31.7 - 28.3)


def test_calculate_delta_ref_defaults_to_zero_per_element(nmr_data):
    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1, "H": 31.7}, delta_ref={"C": 2.0})
    result = calc.calculate(nmr_data)
    assert result["H"]["3"] == pytest.approx(0.0 + 31.7 - 28.3)


# ---------------------------------------------------------------------------
# ChemicalShiftCalculator — filtering
# ---------------------------------------------------------------------------

def test_calculate_skips_elements_without_reference(nmr_data, calc_c_only):
    result = calc_c_only.calculate(nmr_data)
    assert "H" not in result
    assert "N" not in result


def test_calculate_empty_sigma_ref_returns_empty(nmr_data):
    calc = ChemicalShiftCalculator(sigma_ref={})
    assert calc.calculate(nmr_data) == {}


def test_calculate_preserves_atom_indices(nmr_data, calc_c_only):
    result = calc_c_only.calculate(nmr_data)
    assert set(result["C"].keys()) == {"1", "2"}


# ---------------------------------------------------------------------------
# ChemicalShiftCalculator — edge cases
# ---------------------------------------------------------------------------

def test_calculate_negative_shift_when_sigma_mol_greater_than_sigma_ref():
    nmr = NMRData(method="nmr", chemical_shifts={"C": {"1": 200.0}})
    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    result = calc.calculate(nmr)
    assert result["C"]["1"] == pytest.approx(188.1 - 200.0)
    assert result["C"]["1"] < 0


def test_calculate_zero_shift_when_sigma_mol_equals_sigma_ref():
    nmr = NMRData(method="nmr", chemical_shifts={"C": {"1": 188.1}})
    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    result = calc.calculate(nmr)
    assert result["C"]["1"] == pytest.approx(0.0)


def test_calculate_empty_chemical_shifts_returns_empty():
    nmr = NMRData(method="nmr", chemical_shifts={})
    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    assert calc.calculate(nmr) == {}


def test_calculate_negative_sigma_ref():
    nmr = NMRData(method="nmr", chemical_shifts={"Pd": {"1": -200.0}})
    calc = ChemicalShiftCalculator(sigma_ref={"Pd": -277.67})
    result = calc.calculate(nmr)
    assert result["Pd"]["1"] == pytest.approx(-277.67 - (-200.0))


def test_calculate_all_atoms_of_same_element(nmr_data):
    calc = ChemicalShiftCalculator(sigma_ref={"H": 31.7})
    result = calc.calculate(nmr_data)
    assert result["H"]["3"] == pytest.approx(31.7 - 28.3)
    assert result["H"]["4"] == pytest.approx(31.7 - 29.1)
    assert len(result["H"]) == 2


# ---------------------------------------------------------------------------
# WorkflowAnalyzer.apply_chemical_shifts
# ---------------------------------------------------------------------------

def _make_analyzer_with_jobs(tmp_path: Path, jobs: dict) -> WorkflowAnalyzer:
    """Helper that creates a WorkflowAnalyzer and injects pre-built JobResults."""
    analyzer = WorkflowAnalyzer(root_dir=tmp_path)
    analyzer.jobs = jobs
    return analyzer


def test_apply_chemical_shifts_populates_delta_shifts(tmp_path):
    nmr = NMRData(method="nmr", chemical_shifts={"C": {"1": 120.5}})
    job = JobResults(job_id="job_0002", job_type="nmr", nmr=nmr)
    analyzer = _make_analyzer_with_jobs(tmp_path, {"job_0002": job})

    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    analyzer.apply_chemical_shifts(calc)

    assert analyzer.jobs["job_0002"].nmr.delta_shifts["C"]["1"] == pytest.approx(188.1 - 120.5)


def test_apply_chemical_shifts_skips_jobs_without_nmr(tmp_path):
    job = JobResults(job_id="job_0000", job_type="optimization", electronic_energy=-100.5)
    analyzer = _make_analyzer_with_jobs(tmp_path, {"job_0000": job})

    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    analyzer.apply_chemical_shifts(calc)

    assert job.nmr is None


def test_apply_chemical_shifts_skips_unmatched_elements(tmp_path):
    nmr = NMRData(method="nmr", chemical_shifts={"N": {"1": 100.0}})
    job = JobResults(job_id="job_0002", job_type="nmr", nmr=nmr)
    analyzer = _make_analyzer_with_jobs(tmp_path, {"job_0002": job})

    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    analyzer.apply_chemical_shifts(calc)

    assert analyzer.jobs["job_0002"].nmr.delta_shifts == {}


def test_apply_chemical_shifts_multiple_jobs(tmp_path):
    nmr_a = NMRData(method="nmr", chemical_shifts={"C": {"1": 120.5}})
    nmr_b = NMRData(method="nmr", chemical_shifts={"C": {"1": 118.0}})
    job_a = JobResults(job_id="job_0002", job_type="nmr", nmr=nmr_a)
    job_b = JobResults(job_id="job_0002", job_type="nmr", nmr=nmr_b)
    analyzer = _make_analyzer_with_jobs(tmp_path, {"conformer_a": job_a, "conformer_b": job_b})

    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    analyzer.apply_chemical_shifts(calc)

    assert analyzer.jobs["conformer_a"].nmr.delta_shifts["C"]["1"] == pytest.approx(188.1 - 120.5)
    assert analyzer.jobs["conformer_b"].nmr.delta_shifts["C"]["1"] == pytest.approx(188.1 - 118.0)


def test_apply_chemical_shifts_does_not_overwrite_shieldings(tmp_path):
    nmr = NMRData(method="nmr", chemical_shifts={"C": {"1": 120.5}})
    job = JobResults(job_id="job_0002", job_type="nmr", nmr=nmr)
    analyzer = _make_analyzer_with_jobs(tmp_path, {"job_0002": job})

    calc = ChemicalShiftCalculator(sigma_ref={"C": 188.1})
    analyzer.apply_chemical_shifts(calc)

    assert analyzer.jobs["job_0002"].nmr.chemical_shifts["C"]["1"] == pytest.approx(120.5)
