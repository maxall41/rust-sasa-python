use pdbtbx::StrictnessLevel;
use pyo3::prelude::*;
use rust_sasa::{SASALevel, SASAResult};

#[pyfunction]
fn calculate_sasa_at_protein_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>) -> PyResult<f32> {
    let (mut pdb, _errors) = pdbtbx::open(
        pdb_path,
        StrictnessLevel::Medium
    ).unwrap();
    let result = rust_sasa::calculate_sasa(&pdb,probe_radius,n_points,SASALevel::Protein).unwrap();
    return match result {
        SASAResult::Protein(v) => Ok(v),
        _ => panic!("This will never run")
    };
}

#[pyfunction]
fn calculate_sasa_at_residue_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>) -> PyResult<Vec<f32>> {
    let (mut pdb, _errors) = pdbtbx::open(
        pdb_path,
        StrictnessLevel::Medium
    ).unwrap();
    let result = rust_sasa::calculate_sasa(&pdb,probe_radius,n_points,SASALevel::Residue).unwrap();
    return match result {
        SASAResult::Residue(v) => {
            let mut temp : Vec<f32> = vec![];
            for value in v {
                temp.push(value.value)
            }
            Ok(temp)
        },
        _ => panic!("This will never run")
    };
}

#[pyfunction]
fn calculate_sasa_at_atom_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>) -> PyResult<Vec<f32>> {
    let (mut pdb, _errors) = pdbtbx::open(
        pdb_path,
        StrictnessLevel::Medium
    ).unwrap();
    let result = rust_sasa::calculate_sasa(&pdb,probe_radius,n_points,SASALevel::Atom).unwrap();
    return match result {
        SASAResult::Atom(v) => {
            Ok(v)
        },
        _ => panic!("This will never run")
    };
}

#[pyfunction]
fn calculate_sasa_at_chain_level(pdb_path: String,probe_radius: Option<f32>,n_points: Option<usize>) -> PyResult<Vec<f32>> {
    let (mut pdb, _errors) = pdbtbx::open(
        pdb_path,
        StrictnessLevel::Medium
    ).unwrap();
    let result = rust_sasa::calculate_sasa(&pdb,probe_radius,n_points,SASALevel::Chain).unwrap();
    return match result {
        SASAResult::Chain(v) => {
            let mut temp : Vec<f32> = vec![];
            for value in v {
                temp.push(value.value)
            }
            Ok(temp)
        },
        _ => panic!("This will never run")
    };
}

#[pymodule]
fn rust_sasa_python(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(calculate_sasa_at_protein_level, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_sasa_at_chain_level, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_sasa_at_residue_level, m)?)?;
    m.add_function(wrap_pyfunction!(calculate_sasa_at_atom_level, m)?)?;
    Ok(())
}
