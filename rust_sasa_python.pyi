"""Type stubs for rust_sasa_python module."""

from typing import List

class Protein:
    """Protein-level SASA results."""

    total: float
    polar: float
    non_polar: float

class Residue:
    """Residue-level SASA results."""

    chain_id: str
    residue_name: str
    residue_number: int
    sasa: float

class Chain:
    """Chain-level SASA results."""

    chain_id: str
    sasa: float

class SASACalculator:
    """SASA calculator with configurable parameters."""

    def __init__(self, pdb_path: str) -> None:
        """Initialize SASA calculator with PDB file path.

        Args:
            pdb_path: Path to PDB or mmCIF file

        """

    def with_probe_radius(self, radius: float) -> SASACalculator:
        """Set the probe radius for SASA calculation.

        Args:
            radius: Probe sphere radius in Ångströms (default: 1.4)

        Returns:
            Self for method chaining

        """

    def with_n_points(self, points: int) -> SASACalculator:
        """Set the number of points for surface sampling.

        Args:
            points: Number of sampling points (default: 100)

        Returns:
            Self for method chaining

        """

    def calculate_protein(self) -> Protein:
        """Calculate SASA at the protein level.

        Returns:
            Protein: Object containing total, polar, and non-polar SASA values

        """

    def calculate_chain(self) -> List[Chain]:
        """Calculate SASA at the chain level.

        Returns:
            List[Chain]: List of chain-level SASA results

        """

    def calculate_residue(self) -> List[Residue]:
        """Calculate SASA at the residue level.

        Returns:
            List[Residue]: List of residue-level SASA results

        """

    def calculate_atom(self) -> List[float]:
        """Calculate SASA at the atom level.

        Returns:
            List[float]: List of atom-level SASA values

        """

# Convenience functions
def calculate_protein_sasa(pdb_path: str) -> Protein:
    """Calculate protein-level SASA for a PDB file.

    Args:
        pdb_path: Path to PDB or mmCIF file

    Returns:
        Protein: Object containing total, polar, and non-polar SASA values

    """

def calculate_residue_sasa(pdb_path: str) -> List[Residue]:
    """Calculate residue-level SASA for a PDB file.

    Args:
        pdb_path: Path to PDB or mmCIF file

    Returns:
        List[Residue]: List of residue-level SASA results

    """

def calculate_atom_sasa(pdb_path: str) -> List[float]:
    """Calculate atom-level SASA for a PDB file.

    Args:
        pdb_path: Path to PDB or mmCIF file

    Returns:
        List[float]: List of atom-level SASA values

    """

def calculate_chain_sasa(pdb_path: str) -> List[Chain]:
    """Calculate chain-level SASA for a PDB file.

    Args:
        pdb_path: Path to PDB or mmCIF file

    Returns:
        List[Chain]: List of chain-level SASA results

    """

def calculate_sasa_internal(
    atoms_in: list[tuple[tuple[float, float, float], float, int]],
    probe_radius: float,
    n_points: int,
) -> list[float]:
    """Calculate SASA for a set of atoms.

    Args:
        atoms_in: List of atom coordinates, radii, and indices
        probe_radius: Probe radius
        n_points: Number of points on the sphere

    Returns:
        List[float]: List of SASA values

    """
