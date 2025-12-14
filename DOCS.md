# API Reference

### SASACalculator (Builder Interface)

The main interface for SASA calculations using a builder.

```python
# Create calculator instance
calculator = sasa.SASACalculator("protein.pdb")

# Configure parameters (optional, chainable)
calculator = calculator.with_probe_radius(1.2).with_n_points(2000)

# Perform calculations at different levels of granularity
protein_result = calculator.calculate_protein()
chain_results = calculator.calculate_chain()
residue_results = calculator.calculate_residue()
atom_results = calculator.calculate_atom()
```

#### Constructor

```python
SASACalculator(pdb_path: str)
```

**Parameters:**
- `pdb_path` (str): Path to PDB or mmCIF file

#### Configuration Methods

**`with_probe_radius(radius: float) -> SASACalculator`**
- Set the probe sphere radius in Ångströms (default: 1.4)
- Returns modified calculator for method chaining

**`with_n_points(points: int) -> SASACalculator`**
- Set the number of points for surface sampling (default: 100)
- Returns modified calculator for method chaining

**`with_threads(threads: int) -> SASACalculator`**
- Set the number of threads for parallel processing (default: -1 uses all cores)
- Returns modified calculator for method chaining

**`with_include_hydrogens(include_hydrogens: bool) -> SASACalculator`**
- Set whether to include hydrogen atoms in the calculation (default: False)
- Returns modified calculator for method chaining

**`with_include_hetatms(include_hetatms: bool) -> SASACalculator`**
- Set whether to include HETATM records in the calculation (default: False)
- Returns modified calculator for method chaining

**`with_radii_file(path: str) -> SASACalculator`**
- Customize Van der Waals radii by providing a custom radii file (default: ProtOr by Tsai et al.)
- Returns modified calculator for method chaining

**`with_allow_vdw_fallback(allow: bool) -> SASACalculator`**
- Allow fallback to PDBTBX van der Waals radii when custom radius is not found (default: False)
- Returns modified calculator for method chaining

#### Calculation Methods

**`calculate_protein() -> SASAResult`**
- Calculate SASA for the entire protein
- Returns `SASAResult` with `total`, `polar`, and `non_polar` attributes

**`calculate_residue() -> List[ResidueResult]`**
- Calculate SASA for each residue
- Returns list of `ResidueResult` objects

**`calculate_atom() -> List[float]`**
- Calculate SASA for each atom
- Returns list of SASA values

**`calculate_chain() -> List[ChainResult]`**
- Calculate SASA for each chain
- Returns list of `ChainResult` objects

### Result Objects

**`SASAResult`**
- `total`: Total SASA
- `polar`: Total Polar SASA
- `non_polar`: Total Non-polar SASA

**`ResidueResult`**
- `chain_id`: Chain identifier
- `residue_name`: Three-letter residue code
- `residue_number`: Residue sequence number
- `sasa`: SASA value for each residue

**`ChainResult`**
- `chain_id`: Chain identifier
- `sasa`: SASA value for chain

### Convenience Functions

For simple use cases, convenience functions are available:

```python
# Simple one-liners
protein_result = sasa.calculate_protein_sasa("protein.pdb")
residue_results = sasa.calculate_residue_sasa("protein.pdb")
atom_results = sasa.calculate_atom_sasa("protein.pdb")
chain_results = sasa.calculate_chain_sasa("protein.pdb")
```

### Low-level Functions

For advanced use cases where you need to calculate SASA for custom atom sets:

**`calculate_sasa_internal(atoms_in: list[tuple[tuple[float, float, float], float, int]], probe_radius: float, n_points: int, threads: int) -> list[float]`**

Calculate SASA for a set of atoms directly without PDB file parsing.

**Parameters:**
- `atoms_in`: List of tuples, each containing:
  - Atom coordinates as `(x, y, z)` tuple of floats
  - Atom radius as float
  - Atom index as int
- `probe_radius`: Probe sphere radius in Ångströms
- `n_points`: Number of points for surface sampling
- `threads`: Number of threads for parallel processing (-1 uses all cores)

**Returns:**
- List of SASA values (floats) corresponding to each input atom

**Example:**
```python
import rust_sasa_python as sasa

# Define atoms: [(coordinates, radius, index), ...]
atoms = [
    ((0.0, 0.0, 0.0), 1.7, 0),    # Carbon atom
    ((2.0, 0.0, 0.0), 1.2, 1),    # Hydrogen atom
    ((0.0, 2.0, 0.0), 1.55, 2),   # Nitrogen atom
]

# Calculate SASA
sasa_values = sasa.calculate_sasa_internal(atoms, probe_radius=1.4, n_points=100, threads=-1)
print(f"SASA values: {sasa_values}")
```

**`calculate_sasa_internal_at_residue_level(atoms_in: list[tuple[tuple[float, float, float], float, int]], probe_radius: float, n_points: int, threads: int) -> list[ResidueResult]`**

Calculate SASA at the residue level for a set of atoms directly without PDB file parsing. This function groups atoms by residue ID and calculates the total SASA for each residue.

**Parameters:**
- `atoms_in`: List of tuples, each containing:
  - Atom coordinates as `(x, y, z)` tuple of floats
  - Atom radius as float
  - Residue ID as int (used to group atoms into residues)
- `probe_radius`: Probe sphere radius in Ångströms
- `n_points`: Number of points for surface sampling
- `threads`: Number of threads for parallel processing (-1 uses all cores)

**Returns:**
- List of `ResidueResult` objects with:
  - `chain_id`: Always "UNK" (unknown)
  - `residue_name`: Always "UNK" (unknown)
  - `residue_number`: The residue ID from input
  - `sasa`: Total SASA for the residue

**Example:**
```python
import rust_sasa_python as sasa

# Define atoms with residue grouping
# Residue 1: two atoms, Residue 2: one atom
atoms = [
    ((0.0, 0.0, 0.0), 1.7, 1),    # Carbon in residue 1
    ((1.5, 0.0, 0.0), 1.2, 1),    # Hydrogen in residue 1
    ((0.0, 3.0, 0.0), 1.55, 2),   # Nitrogen in residue 2
]

# Calculate residue-level SASA
residues = sasa.calculate_sasa_internal_at_residue_level(atoms, probe_radius=1.4, n_points=100, threads=-1)
for residue in residues:
    print(f"Residue {residue.residue_number}: {residue.sasa:.2f} Ų")
```

## Usage Examples

### Basic Usage

```python
import rust_sasa_python as sasa

# Create calculator
calc = sasa.SASACalculator("protein.pdb")

# Get protein-level results
result = calc.calculate_protein()
print(f"Total SASA: {result.total:.2f} Ų")
print(f"Hydrophobic ratio: {result.non_polar/result.total:.1%}")
```

### Advanced Configuration

```python
# High-precision calculation
calc = (sasa.SASACalculator("protein.pdb")
        .with_probe_radius(1.2)
        .with_n_points(5000))

result = calc.calculate_protein()
print(f"High-precision SASA: {result.total:.3f} Ų")
```

### Multiple Calculations with Same Parameters

```python
# Configure once, calculate multiple levels
calc = (sasa.SASACalculator("protein.pdb")
        .with_probe_radius(1.4)
        .with_n_points(2000))

protein_result = calc.calculate_protein()
residue_results = calc.calculate_residue()
chain_results = calc.calculate_chain()

print(f"Protein total: {protein_result.total:.2f} Ų")
print(f"Number of residues: {len(residue_results)}")
print(f"Number of chains: {len(chain_results)}")
```

### Working with Residue Results

```python
calc = sasa.SASACalculator("protein.pdb")
residues = calc.calculate_residue()

# Find highly exposed residues
exposed = [r for r in residues if r.sasa > 100.0]
print(f"Found {len(exposed)} highly exposed residues")

for residue in exposed[:5]:
    print(f"{residue.chain_id}:{residue.residue_name}{residue.residue_number} = {residue.sasa:.2f} Ų")

# Group by chain
from collections import defaultdict
by_chain = defaultdict(list)
for residue in residues:
    by_chain[residue.chain_id].append(residue.sasa)

for chain_id, sasa_values in by_chain.items():
    print(f"Chain {chain_id}: {sum(sasa_values):.2f} Ų total")
```

### Parameter Comparison

```python
# Compare different probe radii
radii = [1.0, 1.4, 1.8, 2.2]
for radius in radii:
    calc = sasa.SASACalculator("protein.pdb").with_probe_radius(radius)
    result = calc.calculate_protein()
    print(f"Radius {radius}: {result.total:.2f} Ų")
```

### Advanced Configuration Examples

```python
# Include hydrogen atoms in calculation
calc = (sasa.SASACalculator("protein.pdb")
        .with_include_hydrogens(True)
        .with_n_points(2000))
result = calc.calculate_protein()
print(f"SASA with hydrogens: {result.total:.2f} Ų")

# Include HETATM records (ligands, ions, water)
calc = (sasa.SASACalculator("protein.pdb")
        .with_include_hetatms(True))
result = calc.calculate_protein()
print(f"SASA with HETATMs: {result.total:.2f} Ų")

# Control parallelization
calc = (sasa.SASACalculator("protein.pdb")
        .with_threads(4))  # Use 4 threads
result = calc.calculate_protein()

# Use custom Van der Waals radii file
calc = (sasa.SASACalculator("protein.pdb")
        .with_radii_file("custom_radii.txt")
        .with_allow_vdw_fallback(True))  # Fallback to default if not found
result = calc.calculate_protein()
print(f"SASA with custom radii: {result.total:.2f} Ų")
```

### Error Handling

```python
try:
    calc = sasa.SASACalculator("nonexistent.pdb")
    result = calc.calculate_protein()
except RuntimeError as e:
    print(f"Error: {e}")
    # Handle file not found, parsing errors, etc.
```

## Migration from Old API

If you're upgrading from the previous version:

### Old API (deprecated):
```python
# Old way - function with optional parameters
total, non_polar, polar = calculate_sasa_at_protein_level(
    "protein.pdb", probe_radius=1.4, n_points=960
)
```

### New API (recommended):

**Option 1: Builder pattern (recommended for complex usage)**
```python
# New way - builder pattern
calc = (sasa.SASACalculator("protein.pdb")
        .with_probe_radius(1.4)
        .with_n_points(960))
result = calc.calculate_protein()
print(result.total, result.non_polar, result.polar)
```

**Option 2: Convenience functions (for simple usage)**
```python
# New way - simple function
result = sasa.calculate_protein_sasa("protein.pdb")
print(result.total, result.non_polar, result.polar)
```
