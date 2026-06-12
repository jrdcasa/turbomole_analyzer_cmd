from pathlib import Path
from typing import Optional
from turbomole_analyzer.parsers.base import BaseParser
from turbomole_analyzer.models.results import MolecularGeometry, Atom


class CoordParser(BaseParser):
    def parse(self, file_path: Path) -> Optional[MolecularGeometry]:
        if not file_path.exists():
            return None

        atoms = []
        bohr_to_ang = 0.5291772109

        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()

                    # Skip empty lines, turbomole control tags, and comments
                    if not line or line.startswith('$') or line.startswith('#'):
                        continue

                    parts = line.split()
                    if len(parts) >= 4:
                        # Defensive check: if any of the first 3 elements is a comment or non-numeric, skip
                        if parts[0].startswith('#') or parts[1].startswith('#') or parts[2].startswith('#'):
                            continue

                        try:
                            # Safely parse only the first 3 columns as coordinates
                            x_val = float(parts[0]) * bohr_to_ang
                            y_val = float(parts[1]) * bohr_to_ang
                            z_val = float(parts[2]) * bohr_to_ang

                            # Clean the element string (remove numbers or markers like 'c 1' -> 'C')
                            element_str = parts[3].strip().split('#')[0]  # Remove inline comments if any
                            # Keep only alphabetic characters for the element symbol
                            element_clean = ''.join(
                                [char for char in element_str if char.isalpha()]).lower().capitalize()

                            atoms.append(Atom(
                                x=x_val,
                                y=y_val,
                                z=z_val,
                                element=element_clean
                            ))
                        except ValueError:
                            # Individual line parsing failure shouldn't crash the whole file read
                            continue

            if not atoms:
                return None

            return MolecularGeometry(atoms=atoms, num_atoms=len(atoms))

        except Exception as e:
            print(f"Warning: Failed to parse coord file {file_path}. Error: {e}")
            return None
