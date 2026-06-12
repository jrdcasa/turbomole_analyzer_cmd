import re
from pathlib import Path
from typing import List, Dict, Any


class GradientParser:
    """Parses TURBOMOLE 'gradient' files to reconstruct optimization trajectories."""

    def parse_trajectory(self, file_path: Path) -> List[List[Dict[str, Any]]]:
        """
        Parses the gradient file and returns a list of frames.
        Each frame is a list of atoms with elements and coordinates (in Angstroms).
        """
        if not file_path.exists():
            return []

        frames = []
        bohr_to_ang = 0.5291772109

        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()

            current_frame_atoms = []
            reading_coords = False

            for line in lines:
                line = line.strip()

                # Check for cycle initialization or data boundaries
                if "cycle =" in line.lower() or line.startswith("$gradient"):
                    if current_frame_atoms:
                        frames.append(current_frame_atoms)
                        current_frame_atoms = []
                    reading_coords = True
                    continue

                if line.startswith("$end") or line.startswith("$"):
                    if current_frame_atoms:
                        frames.append(current_frame_atoms)
                        current_frame_atoms = []
                    reading_coords = False
                    continue

                if reading_coords and line:
                    parts = line.split()
                    # Turbomole gradient coordinate lines have 4 components or multiples of 3.
                    # Standard geometry line in gradient: X Y Z Element
                    if len(parts) >= 4:
                        # Defensive parsing to ensure we do not capture force/gradient rows
                        # Coordinate rows end with an alphabetical element character symbol
                        element_raw = parts[3].strip().split('#')[0]
                        element_clean = ''.join([c for c in element_raw if c.isalpha()]).lower().capitalize()

                        if element_clean:
                            try:
                                current_frame_atoms.append({
                                    "element": element_clean,
                                    "x": float(parts[0].replace('D', 'E')) * bohr_to_ang,
                                    "y": float(parts[1].replace('D', 'E')) * bohr_to_ang,
                                    "z": float(parts[2].replace('D', 'E')) * bohr_to_ang,
                                })
                            except ValueError:
                                continue

            # Append trailing frame if present
            if current_frame_atoms:
                frames.append(current_frame_atoms)

        except Exception as e:
            print(f"Warning: Failed to parse gradient file {file_path}. Error: {e}")

        return frames