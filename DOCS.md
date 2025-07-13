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
- Set the number of points for surface sampling (default: 1000)
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
