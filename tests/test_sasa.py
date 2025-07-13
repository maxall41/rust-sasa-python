import time
from pathlib import Path

import pytest
import rust_sasa_python as sasa

# Path to the example CIF file
EXAMPLE_CIF = Path(__file__).parent.parent / "example.cif"


def test_protein_calculation():
    """Test protein-level SASA calculation."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF))
    result = calculator.calculate_protein()

    assert isinstance(result, sasa.Protein)
    assert result.total > 0
    assert result.polar >= 0
    assert result.non_polar >= 0
    assert abs(result.total - (result.polar + result.non_polar)) < 0.001


def test_chain_calculation():
    """Test chain-level SASA calculation."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF))
    results = calculator.calculate_chain()

    assert isinstance(results, list)
    assert len(results) > 0

    for chain in results:
        assert isinstance(chain, sasa.Chain)
        assert isinstance(chain.chain_id, str)
        assert chain.sasa >= 0


def test_residue_calculation():
    """Test residue-level SASA calculation."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF))
    results = calculator.calculate_residue()

    assert isinstance(results, list)
    assert len(results) > 0

    for residue in results:
        assert isinstance(residue, sasa.Residue)
        assert isinstance(residue.chain_id, str)
        assert isinstance(residue.residue_name, str)
        assert isinstance(residue.residue_number, int)
        assert residue.sasa >= 0


def test_atom_calculation():
    """Test atom-level SASA calculation."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF))
    results = calculator.calculate_atom()

    assert isinstance(results, list)
    assert len(results) > 0

    for atom_sasa in results:
        assert isinstance(atom_sasa, float)
        assert atom_sasa >= 0


def test_builder_pattern():
    """Test configuration using builder pattern."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF)).with_probe_radius(1.4).with_n_points(500)
    result = calculator.calculate_protein()

    assert result.total > 0
    assert "probe_radius=1.4" in repr(calculator)
    assert "n_points=500" in repr(calculator)


def test_different_configurations():
    """Test that different configurations give different results."""
    calc1 = sasa.SASACalculator(str(EXAMPLE_CIF)).with_probe_radius(1.2)
    calc2 = sasa.SASACalculator(str(EXAMPLE_CIF)).with_probe_radius(1.6)

    result1 = calc1.calculate_protein()
    result2 = calc2.calculate_protein()

    # Different probe radii should give different results
    assert result1.total != result2.total


def test_convenience_functions():
    """Test all convenience functions work."""
    protein = sasa.calculate_protein_sasa(str(EXAMPLE_CIF))
    chains = sasa.calculate_chain_sasa(str(EXAMPLE_CIF))
    residues = sasa.calculate_residue_sasa(str(EXAMPLE_CIF))
    atoms = sasa.calculate_atom_sasa(str(EXAMPLE_CIF))

    assert isinstance(protein, sasa.Protein)
    assert isinstance(chains, list) and len(chains) > 0
    assert isinstance(residues, list) and len(residues) > 0
    assert isinstance(atoms, list) and len(atoms) > 0


def test_consistency_between_granularities():
    """Test that different granularities give consistent totals."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF)).with_n_points(500)

    protein_result = calculator.calculate_protein()
    chain_results = calculator.calculate_chain()
    residue_results = calculator.calculate_residue()
    atom_results = calculator.calculate_atom()

    # Sum of parts should approximately equal whole
    chain_total = sum(chain.sasa for chain in chain_results)
    residue_total = sum(residue.sasa for residue in residue_results)
    atom_total = sum(atom_results)

    tolerance = 0.1
    assert abs(protein_result.total - chain_total) < tolerance
    assert abs(protein_result.total - residue_total) < tolerance
    assert abs(protein_result.total - atom_total) < tolerance


def test_convenience_vs_calculator():
    """Test that convenience functions give same results as calculator."""
    # Using calculator
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF))
    calc_protein = calculator.calculate_protein()

    # Using convenience function
    conv_protein = sasa.calculate_protein_sasa(str(EXAMPLE_CIF))

    # Results should be identical
    assert abs(calc_protein.total - conv_protein.total) < 0.001


def test_nonexistent_file():
    """Test handling of nonexistent file."""
    with pytest.raises(Exception):
        calculator = sasa.SASACalculator("nonexistent_file.pdb")
        calculator.calculate_protein()


def test_calculation_speed():
    """Test that calculations complete quickly."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF)).with_n_points(1000)

    start_time = time.time()
    result = calculator.calculate_protein()
    end_time = time.time()

    calculation_time = end_time - start_time
    assert calculation_time < 2.0  # Should be very fast
    assert result.total > 0


def test_repr_methods():
    """Test that repr methods work for all objects."""
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF))

    # Test calculator repr
    assert "SASACalculator(" in repr(calculator)

    # Test result reprs
    protein = calculator.calculate_protein()
    chains = calculator.calculate_chain()
    residues = calculator.calculate_residue()

    assert "Protein(" in repr(protein)
    if chains:
        assert "Chain(" in repr(chains[0])
    if residues:
        assert "Residue(" in repr(residues[0])


def test_workflow_example():
    """Test a realistic workflow."""
    # Quick calculation
    quick_result = sasa.calculate_protein_sasa(str(EXAMPLE_CIF))
    assert quick_result.total > 0

    # Detailed analysis
    calculator = sasa.SASACalculator(str(EXAMPLE_CIF)).with_probe_radius(1.4).with_n_points(500)

    protein = calculator.calculate_protein()
    chains = calculator.calculate_chain()
    residues = calculator.calculate_residue()

    # All should be valid
    assert protein.total > 0
    assert len(chains) >= 1
    assert len(residues) >= 1

    # Should have expected structure
    if residues:
        first_residue = residues[0]
        assert hasattr(first_residue, "chain_id")
        assert hasattr(first_residue, "residue_name")
        assert hasattr(first_residue, "residue_number")
        assert first_residue.sasa >= 0
