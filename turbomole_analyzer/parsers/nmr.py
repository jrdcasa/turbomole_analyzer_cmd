import re
from pathlib import Path
from typing import Dict, Optional


class NMRParser:
    """Parses TURBOMOLE NMR calculation outputs (mpshift.out)."""

    def parse(self, file_path: Path) -> Optional[Dict[str, Dict[str, float]]]:
        if not file_path.exists():
            return None

        shifts: Dict[str, Dict[str, float]] = {}
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()

            current_element = None
            current_index = None

            for line in lines:
                line_lower = line.lower().strip()

                # Match line: "ATOM  pd   1      ISOTROPIC:      -277.672910       ANISOTROPIC:       338.359117"
                if "atom" in line_lower and "isotropic: " in line_lower:
                    parts = line.split()
                    if len(parts) >= 3:
                        current_index = str(parts[2])  # Capture atom index (e.g., "1")
                        raw_element = parts[1]
                        current_element = ''.join([c for c in raw_element if c.isalpha()]).capitalize()
                        isotropic = float(parts[4])

                        if current_element not in shifts:
                            shifts[current_element] = {}
                        shifts[current_element][current_index] = isotropic

        except Exception as e:
            print(f"Warning: Failed to parse NMR shieldings from {file_path.name}. Error: {e}")

        return shifts if shifts else None