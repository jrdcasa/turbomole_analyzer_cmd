import re
from pathlib import Path
from typing import Optional
from turbomole_analyzer.parsers.base import BaseParser

class FrequencyParser(BaseParser):
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

            # Patterns tried in order; most specific first to avoid false positives.
            # aoforce.out: "zero point vibrational energy :    0.541203 Hartree"
            # aoforce.out: "vibrational zero-point energy:    0.541203 Ha"
            # job.last:    "Zero Point Energy :    0.0750000 Hartree"
            patterns = [
                r"(?:zero[- ]point\s+vibrational\s+energy)\s*:\s*([\-\d\.]+)",
                r"(?:vibrational\s+zero[- ]point\s+energy)\s*:\s*([\-\d\.]+)",
                r"(?:zero[- ]point\s+energy)\s*:\s*([\-\d\.]+)",
            ]

            for pattern in patterns:
                match = re.search(pattern, content, re.IGNORECASE)
                if match:
                    zpe_value = float(match.group(1))
                    break

        except Exception as e:
            print(f"Warning: Failed to parse frequencies from {file_path.name}. Error: {e}")

        return zpe_value