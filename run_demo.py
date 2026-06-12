import tempfile
from pathlib import Path
from turbomole_analyzer.analyzers.workflow import WorkflowAnalyzer

def main():
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)
        trabajo01 = root / "trabajo01"
        
        j0 = trabajo01 / "job_0000"
        j1 = trabajo01 / "job_0001"
        j2 = trabajo01 / "job_0002"
        
        for path in [j0, j1, j2]:
            path.mkdir(parents=True)
            
        with open(j0 / "energy", "w") as f:
            f.write("$energy\n 1  -100.5000000\n$end")
        with open(j0 / "coord", "w") as f:
            f.write("$coord\n 0.0 0.0 0.0 c\n$end")
            
        with open(j1 / "energy", "w") as f:
            f.write("$energy\n 1  -100.4998000\n$end")
        with open(j1 / "job.last", "w") as f:
            f.write("Zero Point Energy : 0.0750000 Hartree\n")
            
        with open(j2 / "energy", "w") as f:
            f.write("$energy\n 1  -100.5015000\n$end")
        with open(j2 / "job.last", "w") as f:
            f.write("Atom   1 C   isotropic shielding =    120.0 ppm\n")

        analyzer = WorkflowAnalyzer(root_dir=root)
        analyzer.process_project("trabajo01")
        
        analyzer.print_individual_reports()
        analyzer.print_relative_energy_report(reference_job_id="job_0000")

if __name__ == "__main__":
    main()
