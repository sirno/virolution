use crate::core::Haplotype;
use crate::readwrite::HaplotypeIO;
use crate::references::HaplotypeRef;
use pyo3::prelude::*;

#[pyclass(name = "Haplotype")]
struct PyHaplotype {
    reference: HaplotypeRef,
}

#[pymethods]
impl PyHaplotype {
    #[new]
    fn new(path: &str) -> Self {
        let wildtype = Haplotype::load_wildtype(path);
        Self {
            reference: wildtype.expect("Failed to load wildtype"),
        }
    }

    fn get_string(&self) -> String {
        self.reference.get_string()
    }

    fn __str__(&self) -> String {
        self.reference.get_string()
    }

    fn create_descendant(
        &self,
        positions: Vec<usize>,
        changes: Vec<u8>,
        generation: usize,
    ) -> Self {
        let optional_changes = changes.into_iter().map(Some).collect();
        let descendant = self
            .reference
            .create_descendant(positions, optional_changes, generation);
        Self {
            reference: descendant,
        }
    }
}

#[pymodule]
fn virolution(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_class::<PyHaplotype>()?;
    Ok(())
}
