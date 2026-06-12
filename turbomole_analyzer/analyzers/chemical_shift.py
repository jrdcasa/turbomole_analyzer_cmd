from typing import Dict, Optional
from turbomole_analyzer.models.results import NMRData


def parse_element_values(arg: str) -> Dict[str, float]:
    """Parse 'C=188.1,H=31.7' into {'C': 188.1, 'H': 31.7}."""
    result: Dict[str, float] = {}
    for item in arg.split(","):
        item = item.strip()
        if not item:
            continue
        element, value = item.split("=", 1)
        result[element.strip().capitalize()] = float(value.strip())
    return result


class ChemicalShiftCalculator:
    """Computes NMR chemical shifts from isotropic shieldings.

    Formula: delta_mol = delta_ref + sigma_ref - sigma_mol
    """

    def __init__(
        self,
        sigma_ref: Dict[str, float],
        delta_ref: Optional[Dict[str, float]] = None,
    ):
        self.sigma_ref = sigma_ref
        self.delta_ref = delta_ref or {}

    def calculate(self, nmr_data: NMRData) -> Dict[str, Dict[str, float]]:
        """Return {element: {atom_idx: delta_mol}} for elements with a defined sigma_ref."""
        result: Dict[str, Dict[str, float]] = {}
        for element, atom_shieldings in nmr_data.chemical_shifts.items():
            if element not in self.sigma_ref:
                continue
            sigma_r = self.sigma_ref[element]
            delta_r = self.delta_ref.get(element, 0.0)
            result[element] = {
                atom_idx: delta_r + sigma_r - sigma_mol
                for atom_idx, sigma_mol in atom_shieldings.items()
            }
        return result
