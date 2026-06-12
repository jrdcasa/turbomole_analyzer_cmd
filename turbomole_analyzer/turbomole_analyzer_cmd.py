import argparse
from pathlib import Path
from typing import List, Dict, Any
import pandas as pd

from turbomole_analyzer.analyzers.workflow import WorkflowAnalyzer
from turbomole_analyzer.parsers.gradient import GradientParser
from turbomole_analyzer.utils.io_helpers import (
    TeeLogger,
    write_optimization_movie,
    write_vmd_tcl_script,
    append_geometry_to_movie_xyz
)


def parse_arguments():
    """
    Configures and parses the command-line interface arguments.
    Enforces mutual exclusivity for input directory routing.
    """
    parser = argparse.ArgumentParser(
        description="Automatic analysis framework for TURBOMOLE calculations."
    )

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-r", "--root",
        type=str,
        help="Root directory containing multiple project folders (e.g., rotamers, conformers...)"
    )
    group.add_argument(
        "-d", "--dirs",
        type=str,
        nargs="+",
        help="Space-separated list of explicit project directory paths."
    )

    parser.add_argument(
        "--ref",
        type=str,
        default=None,
        help="Job ID to use as internal reference within each project. If not given, lowest energy is chosen."
    )

    return parser.parse_args()


def run_analysis_pipeline(args, project_paths: List[Path]):
    """
    Executes the analytical pipeline over all selected paths, generating reports,
    individual optimization movies, and consolidating structural and NMR datasets.
    """
    global_summary_data: Dict[str, Dict[str, Any]] = {}

    # Trackers for multi-conformer NMR cross-matrix generation
    # Structure: { job_id: { element: { atom_idx: { conformer_name: shielding_value } } } }
    nmr_matrix_collector: Dict[str, Dict[str, Dict[str, Dict[str, float]]]] = {}

    # 1. Process each individual project directory
    for proj_path in project_paths:
        print("\n" + "=" * 80)
        print(f" PROCESSING PROJECT DIRECTORY: {proj_path.resolve()}")
        print("=" * 80)

        # Initialize analyzer bound to the parent path framework
        analyzer = WorkflowAnalyzer(root_dir=proj_path.parent)

        try:
            analyzer.process_project(proj_path.name)
        except Exception as e:
            print(f"Error processing {proj_path.name}: {e}")
            continue

        if not analyzer.jobs:
            print(f"[INFO] No valid 'job_XXXX' subdirectories found in {proj_path.name}.")
            continue

        # Phase 1 & 2: Local reports
        analyzer.print_individual_reports()

        chosen_ref = args.ref
        if not chosen_ref:
            valid_jobs = {k: v for k, v in analyzer.jobs.items() if v.electronic_energy is not None}
            if valid_jobs:
                chosen_ref = min(valid_jobs, key=lambda k: valid_jobs[k].electronic_energy)
            else:
                chosen_ref = list(analyzer.jobs.keys())[0]

        analyzer.print_relative_energy_report(reference_job_id=chosen_ref)

        # Collect data points for Phase 3 Global Summary Table
        opt_job = analyzer.jobs.get("job_0000")
        freq_job = analyzer.jobs.get("job_0001")

        global_summary_data[proj_path.name] = {
            "opt_energy": opt_job.electronic_energy if opt_job else None,
            "freq_energy": freq_job.total_energy_with_zpe if freq_job and freq_job.total_energy_with_zpe is not None else (
                freq_job.electronic_energy if freq_job else None)
        }

        # FEATURE: Step-by-step optimization movie generation from 'gradient' file
        gradient_file_path = proj_path / "job_0000" / "gradient"
        if gradient_file_path.exists():
            grad_parser = GradientParser()
            frames = grad_parser.parse_trajectory(gradient_file_path)
            if frames:
                movie_filename = Path(f"{proj_path.name}_optimization.xyz")
                write_optimization_movie(
                    trajectory_frames=frames,
                    output_path=movie_filename,
                    title=proj_path.name
                )
        else:
            print(f"[INFO] No 'gradient' trajectory history file found under {proj_path.name}/job_0000/")

        # =====================================================================
        # FEATURE: Populate cross-conformational data maps for NMR matrices
        # =====================================================================
        for job_id, job_data in analyzer.jobs.items():
            # Check if this job actually contains successfully parsed NMR shifts
            if job_data.nmr and job_data.nmr.chemical_shifts:
                # We group by the real directory name (e.g., 'job_0003') to separate methods
                if job_id not in nmr_matrix_collector:
                    nmr_matrix_collector[job_id] = {}

                for element, atom_dict in job_data.nmr.chemical_shifts.items():
                    if element not in nmr_matrix_collector[job_id]:
                        nmr_matrix_collector[job_id][element] = {}

                    for atom_idx, shielding_val in atom_dict.items():
                        if atom_idx not in nmr_matrix_collector[job_id][element]:
                            nmr_matrix_collector[job_id][element][atom_idx] = {}

                        # Store the shielding mapped to the specific conformer/project folder name
                        nmr_matrix_collector[job_id][element][atom_idx][proj_path.name] = shielding_val

        # =====================================================================
        # FEATURE: Export structured multi-conformer NMR matrices (.dat files)
        # =====================================================================
        conformer_names = [p.name for p in project_paths]

        if not nmr_matrix_collector:
            print("[INFO] No NMR data blocks were detected across any of the processed directories.")

        for job_folder_name, elements_dict in nmr_matrix_collector.items():
            for element, atoms_dict in elements_dict.items():
                # Dynamically names the file based on the actual directory: e.g., H_nmr_job_0003.dat
                output_dat_name = Path(f"{element}_nmr_{job_folder_name}.dat")

                try:
                    with open(output_dat_name, "w", encoding="utf-8") as dat_file:
                        # Write header block compatible with plotting software
                        headers = " ".join([f"shielding_{conf}" for conf in conformer_names])
                        dat_file.write(f"# Atom_label {headers}\n")

                        # Sort atom indices numerically (converting string keys to int for sorting)
                        sorted_atom_indices = sorted(atoms_dict.keys(), key=lambda x: int(x))

                        for atom_idx in sorted_atom_indices:
                            row_values = []
                            for conf in conformer_names:
                                val = atoms_dict[atom_idx].get(conf, "None")
                                if isinstance(val, float):
                                    row_values.append(f"{val:14.4f}")
                                else:
                                    row_values.append(f"{val:>14}")

                            values_str = " ".join(row_values)
                            dat_file.write(f"{atom_idx:<10} {values_str}\n")

                    print(f"[INFO] Exported multi-conformer NMR matrix to: {output_dat_name.resolve()}")
                except Exception as e:
                    print(f"Warning: Failed to write NMR matrix file {output_dat_name.name}. Error: {e}")


    # Phase 3: Global Project Summary across all parsed directories
    if len(global_summary_data) > 1:
        print("\n" + "=" * 95)
        print("           PHASE 3: GLOBAL CROSS-PROJECT SUMMARY (ALL ROTAMERS COMPARED)")
        print("=" * 95)

        valid_opt_energies = [v["opt_energy"] for v in global_summary_data.values() if v["opt_energy"] is not None]
        valid_freq_energies = [v["freq_energy"] for v in global_summary_data.values() if v["freq_energy"] is not None]

        global_min_opt = min(valid_opt_energies) if valid_opt_energies else None
        global_min_freq = min(valid_freq_energies) if valid_freq_energies else None

        hartree_to_kcal = 627.509
        global_rows = []

        for proj_name, data in global_summary_data.items():
            opt_val = data["opt_energy"]
            freq_val = data["freq_energy"]

            delta_e = (opt_val - global_min_opt) * hartree_to_kcal if opt_val and global_min_opt else None
            delta_e_vib = (freq_val - global_min_freq) * hartree_to_kcal if freq_val and global_min_freq else None

            global_rows.append({
                "JOB ID": proj_name,
                "Optimization_calc (Ha)": f"{opt_val:.6f}" if opt_val else "N/A",
                "FEnergy+Vibration (Ha)": f"{freq_val:.6f}" if freq_val else "N/A",
                "Delta E(kcal/mol)": f"{delta_e:.3f}" if delta_e is not None else "N/A",
                "Delta E+Evib(kcal/mol)": f"{delta_e_vib:.3f}" if delta_e_vib is not None else "N/A"
            })

        df_global = pd.DataFrame(global_rows)
        print(df_global.to_string(index=False))
        print("=" * 95)


def main():
    """Main CLI execution block."""
    args = parse_arguments()
    project_paths: List[Path] = []

    # Path discovery logic
    if args.root:
        root_path = Path(args.root)
        if not root_path.exists() or not root_path.is_dir():
            print(f"Error: Root directory '{args.root}' does not exist.")
            return
        for p in sorted(root_path.iterdir()):
            if p.is_dir():
                has_jobs = any(child.is_dir() and child.name.startswith("job_") for child in p.iterdir())
                if has_jobs:
                    project_paths.append(p)
    elif args.dirs:
        for folder in args.dirs:
            p = Path(folder)
            if p.exists() and p.is_dir():
                project_paths.append(p)

    if not project_paths:
        print("No valid project directories containing 'job_*' subfolders were found.")
        return

    # Automatically write the VMD wrapper visualization script
    vmd_script_path = Path("view_trajectory.tcl")
    write_vmd_tcl_script(vmd_script_path)

    # Trigger dual output capture context manager (Stdout + Log File)
    log_file = Path("analyzer_summary.log")
    with TeeLogger(log_file):
        print(f"[LOGGER] Initializing double stdout stream. File output: {log_file.resolve()}")
        run_analysis_pipeline(args, project_paths)


if __name__ == "__main__":
    main()