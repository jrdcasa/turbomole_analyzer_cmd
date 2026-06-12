import sys
from pathlib import Path
from turbomole_analyzer.models.results import MolecularGeometry

class TeeLogger:
    """Custom context manager to redirect stdout to both console and a log file."""
    def __init__(self, log_path: Path):
        self.log_path = log_path
        self.file = None
        self.terminal = sys.stdout

    def __enter__(self):
        self.file = open(self.log_path, "w", encoding="utf-8")
        sys.stdout = self
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self.terminal
        if self.file:
            self.file.close()

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()


def append_geometry_to_movie_xyz(geometry: MolecularGeometry, title: str, output_path: Path):
    """
    Appends a MolecularGeometry block into a single multi-structure (movie) XYZ file.
    """
    try:
        # Open in 'a' (append) mode so frames are concatenated
        with open(output_path, "a", encoding="utf-8") as f:
            # 1. Number of atoms for this frame
            f.write(f"{geometry.num_atoms}\n")
            # 2. Frame comment title (shows up in visualizers as the frame name)
            f.write(f"Frame: {title}\n")
            # 3. Coordinates block
            for atom in geometry.atoms:
                f.write(f"{atom.element:<3} {atom.x:14.8f} {atom.y:14.8f} {atom.z:14.8f}\n")
    except Exception as e:
        print(f"Warning: Failed to append frame to movie XYZ for {title}. Error: {e}")


def write_optimization_movie(trajectory_frames: list, output_path: Path, title: str):
    """
    Writes a full list of optimization frames into a single movie XYZ file.
    """
    if not trajectory_frames:
        return

    try:
        with open(output_path, "w", encoding="utf-8") as f:
            for index, frame in enumerate(trajectory_frames):
                num_atoms = len(frame)
                f.write(f"{num_atoms}\n")
                f.write(f"{title} - Optimization Step {index}\n")
                for atom in frame:
                    f.write(f"{atom['element']:<3} {atom['x']:14.8f} {atom['y']:14.8f} {atom['z']:14.8f}\n")
        print(f"[INFO] Optimization movie generated: {output_path.name} ({len(trajectory_frames)} steps)")
    except Exception as e:
        print(f"Warning: Failed to write optimization movie for {title}. Error: {e}")

def write_vmd_tcl_script(output_path: Path):
    """
    Generates a standard VMD/TCL script optimized for chemical visualization.
    Sets up a CPK representation, white background, and orthographic projection.
    """
    tcl_content = """# =============================================================================
# VMD/TCL Script to automatically load and format optimization trajectories
# Generated automatically by Turbomole Analyzer
# =============================================================================

if { [llength $argv] < 1 } {
    puts "Error: No XYZ file specified."
    puts "Usage: vmd -e view_trajectory.tcl -args <filename.xyz>"
    quit
}

set input_xyz [lindex $argv 0]
puts "[INFO] Loading trajectory file: $input_xyz"

set mol_id [mol new $input_xyz type xyz waitfor all]

# Global Display Settings
color Display Background white
display projection Orthographic
display depthcue off
display axes off

# Graphical Representation Settings (CPK style with Element Colors)
mol delrep 0 $mol_id
mol representation CPK 1.000000 0.300000 12.000000 10.000000
mol selection all
mol coloring Element
mol material Opaque
mol addrep $mol_id

display resetview
puts "[SUCCESS] Trajectory loaded. Use the VMD animation bar to play the optimization movie."
"""
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(tcl_content)
        print(f"[INFO] Created VMD visualization script at: {output_path.name}")
    except Exception as e:
        print(f"Warning: Failed to create VMD script. Error: {e}")