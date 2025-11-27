# rust-sasa-python

A high-performance **Python library** for computing SASA (Solvent Accessible Surface Area) using Rust as the backend. This library provides a clean, Pythonic interface to the fast RustSASA calculation engine.

## Features

- 🦀 **Powered by [RustSASA](https://github.com/maxall41/RustSASA)**: Leverages Rust's performance and safety.
- ⚡️ **Ludicrous Speed**: **63X** faster than Biopython, **5X** faster than Freesasa.
- 🐍 **Pythonic Interface**: Clean, intuitive API.
- 🔧 **Configurable**: Customizable probe radius and sampling points.
- 📁 **PDB and mmCIF SUPPORT**: Supports both PDB and mmCIF files.

## Installation

```bash
pip install rust-sasa-python
```

## Quick Start

```python
import rust_sasa_python as sasa

# Simple calculation - use convenience function
result = sasa.calculate_protein_sasa("protein.pdb")
print(f"Total SASA: {result.total:.2f}")

# Builder pattern for more control
calculator = (sasa.SASACalculator("protein.pdb")
              .with_probe_radius(1.2)
              .with_n_points(2000))
result = calculator.calculate_protein()
print(f"Total SASA: {result.total:.2f}")
print(f"Polar SASA: {result.polar:.2f}")
print(f"Non-polar SASA: {result.non_polar:.2f}")
```

See [DOCS](https://github.com/maxall41/rust-sasa-python/blob/main/DOCS.md) for more information and API reference.

## Contributing

Contributions are welcome! Please feel free to submit issues and pull requests.

## License

This project is licensed under the MIT License.

## Related Projects

- [RustSASA](https://github.com/maxall41/RustSASA) - The core Rust library.
- [DPXRust](https://github.com/maxall41/DPXRust) - Rust library for DPX calculations.
