from pathlib import Path
from typing import Optional
from turbomole_analyzer.parsers.base import BaseParser

class EnergyParser(BaseParser):
    def parse(self, file_path: Path) -> Optional[float]:
        if not file_path.exists():
            return None
        try:
            with open(file_path, 'r') as f:
                lines = f.readlines()
            for line in reversed(lines):
                line = line.strip()
                if line.startswith('$') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    return float(parts[1])
        except Exception as e:
            print(f"Warning: Failed to parse energy from {file_path}. Error: {e}")
            return None
        return None
