use nalgebra::Point3;
use pdbtbx;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;
use rust_sasa::options::SASAOptions;
use rust_sasa::{Atom, SASAProcessor, calculate_sasa_internal as calculate_sasa_internal_internal};
use rust_sasa::{AtomLevel, ChainLevel, ProteinLevel, ResidueLevel}; // As in long long
use std::collections::HashMap;

#[pyclass]
#[derive(Clone)]
pub struct Protein {
    #[pyo3(get)]
    pub total: f32,
    #[pyo3(get)]
    pub polar: f32,
    #[pyo3(get)]
    pub non_polar: f32,
}

#[pymethods]
impl Protein {
    fn __repr__(&self) -> String {
        format!(
            "Protein(total={:.3}, polar={:.3}, non_polar={:.3})",
            self.total, self.polar, self.non_polar
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Residue {
    #[pyo3(get)]
    pub chain_id: String,
    #[pyo3(get)]
    pub residue_name: String,
    #[pyo3(get)]
    pub residue_number: u32,
    #[pyo3(get)]
    pub sasa: f32,
}

#[pymethods]
impl Residue {
    fn __repr__(&self) -> String {
        format!(
            "Residue(chain_id='{}', residue_name='{}', residue_number={}, sasa={:.3})",
            self.chain_id, self.residue_name, self.residue_number, self.sasa
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub struct Chain {
    #[pyo3(get)]
    pub chain_id: String,
    #[pyo3(get)]
    pub sasa: f32,
}

#[pymethods]
impl Chain {
    fn __repr__(&self) -> String {
        format!("Chain(chain_id='{}', sasa={:.3})", self.chain_id, self.sasa)
    }
}

#[pyclass]
#[derive(Clone)]
pub struct SASACalculator {
    pdb_path: String,
    probe_radius: Option<f32>,
    n_points: Option<usize>,
}

#[pymethods]
impl SASACalculator {
    #[new]
    pub fn new(pdb_path: String) -> Self {
        Self {
            pdb_path,
            probe_radius: None,
            n_points: None,
        }
    }

    /// Set the probe radius for SASA calculation
    ///
    /// Args:
    ///     radius: Probe sphere radius in Ångströms (default: 1.4)
    ///
    /// Returns:
    ///     Self for method chaining
    pub fn with_probe_radius(&mut self, radius: f32) -> Self {
        self.probe_radius = Some(radius);
        self.clone()
    }

    /// Set the number of points for surface sampling
    ///
    /// Args:
    ///     points: Number of sampling points (default: 100)
    ///
    /// Returns:
    ///     Self for method chaining
    pub fn with_n_points(&mut self, points: usize) -> Self {
        self.n_points = Some(points);
        self.clone()
    }

    /// Calculate SASA at the protein level
    ///
    /// Returns:
    ///     Protein: Object containing total, polar, and non-polar SASA values
    pub fn calculate_protein(&self) -> PyResult<Protein> {
        let pdb = self.load_pdb()?;
        let mut options = SASAOptions::<ProteinLevel>::new();

        if let Some(radius) = self.probe_radius {
            options = options.with_probe_radius(radius);
        }

        if let Some(points) = self.n_points {
            options = options.with_n_points(points);
        }

        match options.process(&pdb) {
            Ok(result) => Ok(Protein {
                total: result.global_total,
                polar: result.polar_total,
                non_polar: result.non_polar_total,
            }),
            Err(e) => Err(PyRuntimeError::new_err(format!(
                "Failed to calculate protein SASA: {}",
                e
            ))),
        }
    }

    /// Calculate SASA at the chain level
    ///
    /// Returns:
    ///     List[Chain]: List of chain-level SASA results
    pub fn calculate_chain(&self) -> PyResult<Vec<Chain>> {
        let pdb = self.load_pdb()?;
        let mut options = SASAOptions::<ChainLevel>::new();

        if let Some(radius) = self.probe_radius {
            options = options.with_probe_radius(radius);
        }

        if let Some(points) = self.n_points {
            options = options.with_n_points(points);
        }

        match options.process(&pdb) {
            Ok(results) => {
                let mut chain_results = Vec::new();
                for result in results {
                    chain_results.push(Chain {
                        chain_id: result.name,
                        sasa: result.value,
                    });
                }
                Ok(chain_results)
            }
            Err(e) => Err(PyRuntimeError::new_err(format!(
                "Failed to calculate chain SASA: {}",
                e
            ))),
        }
    }

    /// Calculate SASA at the residue level
    ///
    /// Returns:
    ///     List[Residue]: List of residue-level SASA results
    pub fn calculate_residue(&self) -> PyResult<Vec<Residue>> {
        let pdb = self.load_pdb()?;
        let mut options = SASAOptions::<ResidueLevel>::new();

        if let Some(radius) = self.probe_radius {
            options = options.with_probe_radius(radius);
        }

        if let Some(points) = self.n_points {
            options = options.with_n_points(points);
        }

        match options.process(&pdb) {
            Ok(results) => {
                let mut residue_results = Vec::new();
                for result in results {
                    residue_results.push(Residue {
                        chain_id: result.chain_id,
                        residue_name: result.name,
                        residue_number: result.serial_number as u32,
                        sasa: result.value,
                    });
                }
                Ok(residue_results)
            }
            Err(e) => Err(PyRuntimeError::new_err(format!(
                "Failed to calculate residue SASA: {}",
                e
            ))),
        }
    }

    /// Calculate SASA at the atom level
    ///
    /// Returns:
    ///     List[float]: List of atom-level SASA values
    pub fn calculate_atom(&self) -> PyResult<Vec<f32>> {
        let pdb = self.load_pdb()?;
        let mut options = SASAOptions::<AtomLevel>::new();

        if let Some(radius) = self.probe_radius {
            options = options.with_probe_radius(radius);
        }

        if let Some(points) = self.n_points {
            options = options.with_n_points(points);
        }

        match options.process(&pdb) {
            Ok(results) => Ok(results),
            Err(e) => Err(PyRuntimeError::new_err(format!(
                "Failed to calculate atom SASA: {}",
                e
            ))),
        }
    }

    fn __repr__(&self) -> String {
        let mut parts = vec![format!("SASACalculator(pdb_path='{}')", self.pdb_path)];

        if let Some(radius) = self.probe_radius {
            parts.push(format!("probe_radius={}", radius));
        }

        if let Some(points) = self.n_points {
            parts.push(format!("n_points={}", points));
        }

        if parts.len() > 1 {
            format!(
                "SASACalculator(pdb_path='{}', {})",
                self.pdb_path,
                parts[1..].join(", ")
            )
        } else {
            parts[0].clone()
        }
    }
}

impl SASACalculator {
    fn load_pdb(&self) -> PyResult<pdbtbx::PDB> {
        match pdbtbx::open(&self.pdb_path) {
            Ok((pdb, _errors)) => Ok(pdb),
            Err(e) => Err(PyRuntimeError::new_err(format!(
                "Failed to load PDB file '{}': {:?}",
                self.pdb_path, e
            ))),
        }
    }
}

// Convenience functions for backward compatibility and simple use cases
#[pyfunction]
#[pyo3(signature = (pdb_path))]
pub fn calculate_protein_sasa(pdb_path: String) -> PyResult<Protein> {
    SASACalculator::new(pdb_path).calculate_protein()
}

#[pyfunction]
#[pyo3(signature = (pdb_path))]
pub fn calculate_residue_sasa(pdb_path: String) -> PyResult<Vec<Residue>> {
    SASACalculator::new(pdb_path).calculate_residue()
}

#[pyfunction]
#[pyo3(signature = (pdb_path))]
pub fn calculate_atom_sasa(pdb_path: String) -> PyResult<Vec<f32>> {
    SASACalculator::new(pdb_path).calculate_atom()
}

#[pyfunction]
#[pyo3(signature = (pdb_path))]
pub fn calculate_chain_sasa(pdb_path: String) -> PyResult<Vec<Chain>> {
    SASACalculator::new(pdb_path).calculate_chain()
}

#[pyfunction]
pub fn calculate_sasa_internal(
    atoms_in: Vec<((f32, f32, f32), f32, usize)>,
    probe_radius: f32,
    n_points: usize,
) -> PyResult<Vec<f32>> {
    let atoms: Vec<Atom> = atoms_in
        .into_iter()
        .map(|(pos, radius, index)| Atom {
            position: Point3::new(pos.0, pos.1, pos.2),
            id: index,
            parent_id: None,
            radius,
        })
        .collect();
    Ok(calculate_sasa_internal_internal(
        atoms.as_slice(),
        probe_radius,
        n_points,
        true,
    ))
}

#[pyfunction]
fn calculate_sasa_internal_at_residue_level(
    atoms_in: Vec<((f32, f32, f32), f32, usize)>,
    probe_radius: f32,
    n_points: usize,
) -> PyResult<Vec<Residue>> {
    let atom_sasa = calculate_sasa_internal(atoms_in.clone(), probe_radius, n_points).unwrap();

    let mut residue_groups: HashMap<usize, Vec<f32>> = HashMap::new();

    for (i, atom) in atoms_in.iter().enumerate() {
        let parent_id = atom.2;
        residue_groups
            .entry(parent_id)
            .or_insert_with(Vec::new)
            .push(atom_sasa[i]);
    }

    let mut residue_sasa = Vec::new();
    for (parent_id, sasa_values) in residue_groups {
        let local_sum = sasa_values.iter().sum::<f32>();

        residue_sasa.push(Residue {
            chain_id: "UNK".to_string(),
            residue_name: "UNK".to_string(),
            residue_number: parent_id as u32,
            sasa: local_sum / sasa_values.len() as f32,
        });
    }

    Ok(residue_sasa)
}

#[pymodule]
fn rust_sasa_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<SASACalculator>()?;
    m.add_class::<Protein>()?;
    m.add_class::<Residue>()?;
    m.add_class::<Chain>()?;

    // Convenience functions
    m.add_function(wrap_pyfunction!(calculate_protein_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_residue_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_atom_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_chain_sasa, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_sasa_internal, m)?)?;
    m.add_function(wrap_pyfunction!(
        calculate_sasa_internal_at_residue_level,
        m
    )?)?;

    Ok(())
}
