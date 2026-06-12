import re
from pathlib import Path
from typing import Optional

class FrequencyParser:
    """Parses TURBOMOLE frequency calculation outputs (aoforce.out / job.last)."""

    def parse(self, file_path: Path) -> Optional[float]:
        """
        Extracts the Zero Point Vibrational Energy (ZPE) in Hartree.
        Works with both aoforce.out and job.last files.
        """
        if not file_path.exists():
            return None

        zpe_value = None
        try:
            with open(file_path, 'r') as f:
                content = f.read()

            # Pattern for standard aoforce.out summary blocks
            # Example: "zero point vibrational energy :    0.541203 Hartree"
            # Example: "Zero-point vibrational energy:    0.541203 Ha"
            patterns = [
                r"(?:zero[- ]point\s+vibrational\s+energy)\s*:\s*([\-\d\.]+)",
                r"(?:vibrational\s+zero[- ]point\s+energy)\s*:\s*([\-\d\.]+)"
            ]

            for pattern in patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    zpe_value = float(match.group(1))
                    break

        except Exception as e:
            print(f"Warning: Failed to parse frequencies from {file_path.name}. Error: {e}")

        return zpe_value