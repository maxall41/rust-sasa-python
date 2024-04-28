# rust-sasa-python

rust-sasa-python is a Python library for computing SASA (Solvent Accessible Surface Area) far faster than biopython and other packages.
See [RustSasa](https://github.com/maxall41/RustSASA) for more info.

## Features:
- 🦀 Backend written in Pure Rust
- ⚡️ 3X Faster than Biopython and ~120% faster than Freesasa
- 🐍 Quick drop-in replacement for Biopython

## Usage

You can now utilize RustSasa within Python to speed up your scripts! Take a look at [rust-sasa-python](https://github.com/maxall41/rust-sasa-python)!

Installation:
```
pip install rust-sasa-python
```
Example:
```python
from rust_sasa_python import calculate_sasa_at_residue_level
residue_sasa_values = calculate_sasa_at_residue_level("path_to_pdb_file.pdb") # Also supports mmCIF files!
```
See full docs [here](https://github.com/maxall41/rust-sasa-python/blob/main/DOCS.md)

## Benchmarking
Benchmarks were performed on an M2 Apple Mac with 8GB of RAM and 8 Cores with the protein AF-A0A2K5XT84-F1 (AlphaFold).

- Biopython: ~150ms

- Freesasa: ~90ms

- RustSASA: ~40ms
