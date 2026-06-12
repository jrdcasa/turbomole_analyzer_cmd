from pathlib import Path
from typing import TYPE_CHECKING, Dict
import pandas as pd

from turbomole_analyzer.parsers.coord import CoordParser
from turbomole_analyzer.parsers.energy import EnergyParser
from turbomole_analyzer.parsers.frequency import FrequencyParser
from turbomole_analyzer.parsers.nmr import NMRParser
from turbomole_analyzer.models.results import JobResults, NMRData

if TYPE_CHECKING:
    from turbomole_analyzer.analyzers.chemical_shift import ChemicalShiftCalculator


class WorkflowAnalyzer:
    def __init__(self, root_dir: Path):
        self.root_dir = Path(root_dir)
        self.jobs: Dict[str, JobResults] = {}

    def _detect_job_type(self, job_dir: Path) -> str:
        """
        Heuristically deduces the TURBOMOLE job type based on file presence and contents.
        """
        if (job_dir / "gradient").exists():
            return "optimization"

        if (job_dir / "aoforce.out").exists() and (job_dir / "vibspectrum").exists():
            return "frequency"

        if (job_dir / "mpshift.out").exists() and (job_dir / "shielding").exists():
            return "nmr"

        return "single_point"

    def analyze_job(self, job_dir: Path) -> JobResults:
        """Parses a job directory, looking safely into dedicated module outputs."""
        job_id = job_dir.name
        job_type = self._detect_job_type(job_dir)

        e_elec = EnergyParser().parse(job_dir / "energy")
        geom = CoordParser().parse(job_dir / "coord")

        # 1. CASCADE PARSING FOR FREQUENCIES (aoforce.out -> job.last)
        zpe = None
        if job_type == "frequency":
            aoforce_path = job_dir / "aoforce.out"
            job_last_path = job_dir / "job.last"

            if aoforce_path.exists():
                zpe = FrequencyParser().parse(aoforce_path)
            elif job_last_path.exists():
                zpe = FrequencyParser().parse(job_last_path)

        # 2. CASCADE PARSING FOR NMR (mpshift.out -> job.last)
        nmr_raw = None
        if job_type == "nmr":
            mpshift_path = job_dir / "mpshift.out"
            job_last_path = job_dir / "job.last"

            if mpshift_path.exists():
                nmr_raw = NMRParser().parse(mpshift_path)
            elif job_last_path.exists():
                nmr_raw = NMRParser().parse(job_last_path)

        nmr_data = NMRData(method=job_type, chemical_shifts=nmr_raw) if nmr_raw else None
        e_total_zpe = (e_elec + zpe) if (e_elec is not None and zpe is not None) else None

        # Build results object
        res = JobResults(
            job_id=job_id,
            job_type=job_type,
            electronic_energy=e_elec,
            zpe_correction=zpe,
            total_energy_with_zpe=e_total_zpe,
            geometry=geom,
            nmr=nmr_data
        )

        return res

    def process_project(self, project_dir_name: str):
        """Finds and processes any sorted subdirectories matching 'job_*' pattern."""
        project_path = self.root_dir / project_dir_name
        self.jobs.clear()

        if not project_path.exists():
            raise FileNotFoundError(f"Path {project_path} not found.")

        for job_dir in sorted(project_path.glob("job_*")):
            if job_dir.is_dir():
                job_id = job_dir.name
                result = self.analyze_job(job_dir)

                # Smart energy fallback: If job_0001 lacks an 'energy' file,
                # it inherits the ground-state electronic energy from job_0000
                if result.electronic_energy is None and job_id == "job_0001" and "job_0000" in self.jobs:
                    inherited_energy = self.jobs["job_0000"].electronic_energy
                    zpe = result.zpe_correction
                    e_total = (inherited_energy + zpe) if (inherited_energy is not None and zpe is not None) else None
                    result = JobResults(
                        job_id=result.job_id,
                        job_type=result.job_type,
                        electronic_energy=inherited_energy,
                        zpe_correction=zpe,
                        total_energy_with_zpe=e_total,
                        geometry=result.geometry,
                        nmr=result.nmr
                    )

                self.jobs[job_id] = result

    def apply_chemical_shifts(self, calculator: "ChemicalShiftCalculator") -> None:
        """Compute chemical shifts for all jobs that have NMR shielding data."""
        for job in self.jobs.values():
            if job.nmr and job.nmr.chemical_shifts:
                job.nmr.delta_shifts = calculator.calculate(job.nmr)

    def print_individual_reports(self):
        """Phase 1: Prints a comprehensive structured log per each Job."""
        print("=" * 60)
        print("           PHASE 1: INDIVIDUAL JOB REPORTS")
        print("=" * 60)
        for job_id, job in self.jobs.items():
            print(f"\n[JOB REPORT: {job_id}] ({job.job_type.upper()})")
            print(f"  - Electronic Energy (E_elec) : {job.electronic_energy if job.electronic_energy is not None else 'N/A'} Ha")
            if job.zpe_correction is not None:
                print(f"  - ZPE Correction             : {job.zpe_correction} Ha")
                print(f"  - Total Energy (E_elec+ZPE)  : {job.total_energy_with_zpe} Ha")
            if job.geometry:
                print(f"  - Molecular Geometry        : {job.geometry.num_atoms} atoms detected.")
            if job.nmr and job.nmr.chemical_shifts:
                print(f"  - NMR Shieldings Extracted   : {list(job.nmr.chemical_shifts.keys())}")
            if job.nmr and job.nmr.delta_shifts:
                print(f"  - Chemical Shifts (ppm):")
                for element, atoms in sorted(job.nmr.delta_shifts.items()):
                    for atom_idx, delta in sorted(atoms.items(), key=lambda x: int(x[0])):
                        print(f"      {element}{atom_idx}: {delta:8.3f} ppm")
            print("-" * 40)

    def print_relative_energy_report(self, reference_job_id: str):
        """Phase 2: Compiles the comparative energy table based on a user-defined job reference."""
        print("\n" + "=" * 80)
        print(f"    PHASE 2: RELATIVE ENERGY REPORT (Ref: {reference_job_id})")
        print("=" * 80)

        if reference_job_id not in self.jobs:
            print(f"Error: Reference job '{reference_job_id}' not found in the project.")
            return

        ref_job = self.jobs[reference_job_id]
        ref_e_elec = ref_job.electronic_energy

        # Look for ZPE in the project to define the Freq+ZPE reference
        ref_zpe = None
        for j in self.jobs.values():
            if j.zpe_correction is not None:
                ref_zpe = j.zpe_correction
                break

        ref_total_with_zpe = (ref_e_elec + ref_zpe) if (ref_e_elec is not None and ref_zpe is not None) else None
        hartree_to_kcal = 627.509

        print(f"Reference Electronic Energy (E_el): {ref_e_elec} Ha")
        print(f"Reference ZPE Correction Found    : {ref_zpe if ref_zpe is not None else 'None'} Ha")
        print(f"Reference Total (E_el + ZPE)      : {ref_total_with_zpe if ref_total_with_zpe is not None else 'N/A'} Ha")
        print("-" * 80)
        print(f"{'Job ID':<10} | {'Job Type':<15} | {'Abs Energy (Ha)':<18} | {'dE(elec+ZPE) (kcal/mol)':<22}")
        print("-" * 80)

        for job_id, job in self.jobs.items():
            # 1. Absolute Energy formatting (E_elec in Hartrees/Ha)
            abs_energy_str = f"{job.electronic_energy:16.6f}" if job.electronic_energy is not None else "N/A"

            # 2. Delta Elec + ZPE calculation
            d_zpe_str = "N/A"
            current_zpe = job.zpe_correction if job.zpe_correction is not None else ref_zpe

            if job.electronic_energy is not None and current_zpe is not None and ref_total_with_zpe is not None:
                current_total = job.electronic_energy + current_zpe
                d_zpe = (current_total - ref_total_with_zpe) * hartree_to_kcal
                d_zpe_str = f"{d_zpe:20.3f}"

            print(f"{job_id:<10} | {job.job_type:<15} | {abs_energy_str:<18} | {d_zpe_str:<22}")
        print("=" * 80)