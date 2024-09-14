# Documentation

## Functions

### `calculate_sasa_at_protein_level`

Signature: calculate_sasa_at_protein_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>)

Returns: tuple (GLOBAL_SASA,Non-polar total SASA,Polar total SASA) as f32
Example:

```python
import rust_sasa_python
rust_sasa_python.calculate_sasa_at_protein_level("path_to_pdb_file.pdb") # Also supports mmCIF files!
```

### `calculate_sasa_at_residue_level`

Signature: calculate_sasa_at_residue_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>)

Returns: list of tuples <({CHAIN ID}\_{RESIDUE NAME}_{RESIDUE INDEX},SASA VALUE as f32)>
Example:

```python
import rust_sasa_python
rust_sasa_python.calculate_sasa_at_residue_level("path_to_pdb_file.pdb") # Also supports mmCIF files!
```

### `calculate_sasa_at_atom_level`

Signature: calculate_sasa_at_atom_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>)

Returns: array of SASA values for each atom
Example:

```python
import rust_sasa_python
rust_sasa_python.calculate_sasa_at_atom_level("path_to_pdb_file.pdb") # Also supports mmCIF files!
```

### `calculate_sasa_at_chain_level`

Signature: calculate_sasa_at_chain_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>)

Returns: array of SASA values for each chain
Example:

```python
import rust_sasa_python
rust_sasa_python.calculate_sasa_at_chain_level("path_to_pdb_file.pdb") # Also supports mmCIF files!
```
