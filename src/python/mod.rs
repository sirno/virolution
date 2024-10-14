use crate::core::FitnessProvider;
use crate::core::Haplotype;
use crate::readwrite::HaplotypeIO;
use crate::references::HaplotypeRef;
use pyo3::prelude::*;

#[pyclass(name = "Fitness")]
struct PyFitness {
    fitness_provider: FitnessProvider,
}

#[pymethods]
impl PyFitness {
    #[new]
    #[pyo3(signature = (id, seq, n_symbols, fitness_model))]
    fn create(id: usize, seq: Vec<u8>, n_symbols: usize, fitness_model: &str) -> PyResult<Self> {
        // TODO: Define how to load and interact with fitness from python
        // This should feel as if it was python first...
        unimplemented!()
    }

    fn get_fitness(&self, haplotype: &PyHaplotype) -> f64 {
        self.fitness_provider.get_fitness(&haplotype.reference)
    }
}

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
