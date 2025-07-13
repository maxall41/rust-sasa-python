import rust_sasa_python as sasa

# Simple calculation - use convenience function
calculator = sasa.SASACalculator("example.cif").with_n_points(960)
residues = calculator.calculate_residue()
print(residues)
