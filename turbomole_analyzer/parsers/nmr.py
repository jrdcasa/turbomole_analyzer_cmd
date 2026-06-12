import re
from pathlib import Path
from typing import Dict, Optional
from turbomole_analyzer.parsers.base import BaseParser


class NMRParser(BaseParser):
    """Parses TURBOMOLE NMR calculation outputs (mpshift.out and job.last)."""

    # mpshift.out: "ATOM  pd   1      ISOTROPIC:      -277.672910       ANISOTROPIC: ..."
    _PATTERN_MPSHIFT = re.compile(
        r"atom\s+([a-zA-Z]+)\s+(\d+)\s+isotropic:\s*([\-\d\.]+)",
        re.IGNORECASE,
    )
    # job.last:    "Atom   1 C   isotropic shielding =    120.0 ppm"
    _PATTERN_JOB_LAST = re.compile(
        r"atom\s+(\d+)\s+([a-zA-Z]+)\s+isotropic\s+shielding\s*=\s*([\-\d\.]+)",
        re.IGNORECASE,
    )

    def parse(self, file_path: Path) -> Optional[Dict[str, Dict[str, float]]]:
        if not file_path.exists():
            return None

        shifts: Dict[str, Dict[str, float]] = {}
        try:
            with open(file_path, "r") as f:
                for line in f:
                    m = self._PATTERN_MPSHIFT.search(line)
                    if m:
                        element = "".join(c for c in m.group(1) if c.isalpha()).capitalize()
                        shifts.setdefault(element, {})[m.group(2)] = float(m.group(3))
                        continue

                    m = self._PATTERN_JOB_LAST.search(line)
                    if m:
                        element = "".join(c for c in m.group(2) if c.isalpha()).capitalize()
                        shifts.setdefault(element, {})[m.group(1)] = float(m.group(3))

        except Exception as e:
            print(f"Warning: Failed to parse NMR shieldings from {file_path.name}. Error: {e}")

        return shifts if shifts else None