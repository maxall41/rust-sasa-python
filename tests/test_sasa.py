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
    assert abs(result.total - (result.polar + result.non_polar)) < 1.0


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


def test_with_threads():
    """Test thread configuration."""
    # Test with different thread counts
    calc_single = sasa.SASACalculator(str(EXAMPLE_CIF)).with_threads(1)
    calc_multi = sasa.SASACalculator(str(EXAMPLE_CIF)).with_threads(4)
    calc_auto = sasa.SASACalculator(str(EXAMPLE_CIF)).with_threads(-1)

    result_single = calc_single.calculate_protein()
    result_multi = calc_multi.calculate_protein()
    result_auto = calc_auto.calculate_protein()

    # Results should be identical regardless of thread count
    assert abs(result_single.total - result_multi.total) < 0.01
    assert abs(result_single.total - result_auto.total) < 0.01


def test_with_include_hydrogens():
    """Test hydrogen inclusion configuration."""
    calc_no_h = sasa.SASACalculator(str(EXAMPLE_CIF)).with_include_hydrogens(False)
    calc_with_h = sasa.SASACalculator(str(EXAMPLE_CIF)).with_include_hydrogens(True)

    result_no_h = calc_no_h.calculate_protein()
    result_with_h = calc_with_h.calculate_protein()

    # Both should return valid results
    assert result_no_h.total > 0
    assert result_with_h.total > 0

    # Atom counts should differ if hydrogens are present
    atoms_no_h = calc_no_h.calculate_atom()
    atoms_with_h = calc_with_h.calculate_atom()

    assert len(atoms_no_h) > 0
    assert len(atoms_with_h) > 0


def test_with_include_hetatms():
    """Test HETATM inclusion configuration."""
    calc_no_het = sasa.SASACalculator(str(EXAMPLE_CIF)).with_include_hetatms(False)
    calc_with_het = sasa.SASACalculator(str(EXAMPLE_CIF)).with_include_hetatms(True)

    result_no_het = calc_no_het.calculate_protein()
    result_with_het = calc_with_het.calculate_protein()

    # Both should return valid results
    assert result_no_het.total > 0
    assert result_with_het.total > 0

    # Atom counts may differ if HETATMs are present
    atoms_no_het = calc_no_het.calculate_atom()
    atoms_with_het = calc_with_het.calculate_atom()

    assert len(atoms_no_het) > 0
    assert len(atoms_with_het) > 0


def test_with_allow_vdw_fallback():
    """Test VDW fallback configuration."""
    # Test that fallback can be enabled/disabled
    calc_no_fallback = sasa.SASACalculator(str(EXAMPLE_CIF)).with_allow_vdw_fallback(False)
    calc_with_fallback = sasa.SASACalculator(str(EXAMPLE_CIF)).with_allow_vdw_fallback(True)

    result_no_fallback = calc_no_fallback.calculate_protein()
    result_with_fallback = calc_with_fallback.calculate_protein()

    # Both should return valid results for standard PDB files
    assert result_no_fallback.total > 0
    assert result_with_fallback.total > 0


def test_multiple_config_chaining():
    """Test chaining multiple configuration methods."""
    calculator = (
        sasa.SASACalculator(str(EXAMPLE_CIF))
        .with_probe_radius(1.4)
        .with_n_points(1000)
        .with_threads(2)
        .with_include_hydrogens(False)
        .with_include_hetatms(False)
        .with_allow_vdw_fallback(True)
    )

    # Should be able to calculate at all levels
    protein = calculator.calculate_protein()
    chains = calculator.calculate_chain()
    residues = calculator.calculate_residue()
    atoms = calculator.calculate_atom()

    assert protein.total > 0
    assert len(chains) > 0
    assert len(residues) > 0
    assert len(atoms) > 0


def test_probe_radius_variations():
    """Test various probe radius values."""
    radii = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    results = []

    for radius in radii:
        calc = sasa.SASACalculator(str(EXAMPLE_CIF)).with_probe_radius(radius)
        result = calc.calculate_protein()
        results.append(result.total)

    # All results should be positive
    assert all(r > 0 for r in results)
    # Larger probe radii should generally give smaller SASA values
    # (larger probe can't access smaller crevices)
    assert results[-1] < results[0]


def test_n_points_variations():
    """Test various n_points values for accuracy."""
    points = [50, 100, 500, 1000]
    results = []

    for n in points:
        calc = sasa.SASACalculator(str(EXAMPLE_CIF)).with_n_points(n)
        result = calc.calculate_protein()
        results.append(result.total)

    # Higher point counts should give more consistent results
    assert all(r > 0 for r in results)
    # Results should converge (difference between high point counts should be small)
    assert abs(results[-1] - results[-2]) < abs(results[1] - results[0])


def test_calculate_sasa_internal_with_threads():
    """Test low-level calculate_sasa_internal function with threads parameter."""
    # Define simple atom set
    atoms = [
        ((0.0, 0.0, 0.0), 1.7, 0),
        ((2.0, 0.0, 0.0), 1.2, 1),
        ((0.0, 2.0, 0.0), 1.55, 2),
    ]

    # Test with different thread counts
    sasa_single = sasa.calculate_sasa_internal(atoms, 1.4, 100, 1)
    sasa_multi = sasa.calculate_sasa_internal(atoms, 1.4, 100, 4)
    sasa_auto = sasa.calculate_sasa_internal(atoms, 1.4, 100, -1)

    # Should return same number of values
    assert len(sasa_single) == len(atoms)
    assert len(sasa_multi) == len(atoms)
    assert len(sasa_auto) == len(atoms)

    # Results should be nearly identical
    for i in range(len(atoms)):
        assert abs(sasa_single[i] - sasa_multi[i]) < 0.01
        assert abs(sasa_single[i] - sasa_auto[i]) < 0.01


def test_calculate_sasa_internal_at_residue_level_with_threads():
    """Test low-level residue calculation with threads parameter."""
    # Define atoms grouped by residue
    atoms = [
        ((0.0, 0.0, 0.0), 1.7, 1),
        ((1.5, 0.0, 0.0), 1.2, 1),
        ((0.0, 3.0, 0.0), 1.55, 2),
    ]

    # Test with different thread counts
    residues_single = sasa.calculate_sasa_internal_at_residue_level(atoms, 1.4, 100, 1)
    residues_multi = sasa.calculate_sasa_internal_at_residue_level(atoms, 1.4, 100, 4)
    residues_auto = sasa.calculate_sasa_internal_at_residue_level(atoms, 1.4, 100, -1)

    # Should return same number of residues
    assert len(residues_single) == 2
    assert len(residues_multi) == 2
    assert len(residues_auto) == 2

    # Sort by residue_number to ensure consistent ordering
    residues_single = sorted(residues_single, key=lambda r: r.residue_number)
    residues_multi = sorted(residues_multi, key=lambda r: r.residue_number)
    residues_auto = sorted(residues_auto, key=lambda r: r.residue_number)

    # Results should be nearly identical when matched by residue_number
    for r1, r2, r3 in zip(residues_single, residues_multi, residues_auto):
        assert r1.residue_number == r2.residue_number == r3.residue_number
        assert abs(r1.sasa - r2.sasa) < 0.01
        assert abs(r1.sasa - r3.sasa) < 0.01
