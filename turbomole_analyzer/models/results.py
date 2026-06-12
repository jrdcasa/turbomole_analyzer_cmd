from pydantic import BaseModel, Field
from typing import List, Optional, Dict
import numpy as np

class Atom(BaseModel):
    element: str
    x: float
    y: float
    z: float

class MolecularGeometry(BaseModel):
    atoms: List[Atom] = Field(default_factory=list)

    @property
    def num_atoms(self) -> int:
        return len(self.atoms)

    def get_distance(self, idx1: int, idx2: int) -> float:
        """Calculate the Euclidean distance between two atoms (0-indexed)."""
        if idx1 >= len(self.atoms) or idx2 >= len(self.atoms):
            raise IndexError("Atom index out of range.")
        a1 = self.atoms[idx1]
        a2 = self.atoms[idx2]
        p1 = np.array([a1.x, a1.y, a1.z])
        p2 = np.array([a2.x, a2.y, a2.z])
        return float(np.linalg.norm(p1 - p2))

class NMRData(BaseModel):
    method: str
    # Key: Element (e.g., "H"), Value: Dict mapping atom index string to shielding value
    # Example: {"H": {"1": 31.45, "2": 30.12}, "C": {"3": 135.2}}
    chemical_shifts: Dict[str, Dict[str, float]] = Field(default_factory=dict)
    # Computed chemical shifts: delta_mol = delta_ref + sigma_ref - sigma_mol
    delta_shifts: Dict[str, Dict[str, float]] = Field(default_factory=dict)

class JobResults(BaseModel):
    job_id: str
    job_type: str
    electronic_energy: Optional[float] = None
    zpe_correction: Optional[float] = None
    total_energy_with_zpe: Optional[float] = None
    geometry: Optional[MolecularGeometry] = None
    nmr: Optional[NMRData] = None
