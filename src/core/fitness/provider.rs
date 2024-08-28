use std::fs;
use std::io;
use std::path::Path;

use super::epistasis::EpistasisTable;
use super::init::{FitnessDistribution, FitnessModel};
use super::table::FitnessTable;
use super::utility::UtilityFunction;

use crate::core::haplotype::Symbol;
use crate::errors::VirolutionError;
use crate::references::HaplotypeRef;

#[derive(Clone, Debug)]
pub struct FitnessProvider {
    pub id: usize,
    function: FitnessFunction,
    utility: UtilityFunction,
}

#[derive(Clone, Debug)]
pub enum FitnessFunction {
    NonEpistatic(FitnessTable),
    SimpleEpistatic(FitnessTable, EpistasisTable),
}

impl FitnessProvider {
    pub fn from_model(
        id: usize,
        sequence: &[Symbol],
        n_symbols: usize,
        model: &FitnessModel,
    ) -> Result<Self, VirolutionError> {
        let function = match model.distribution {
            FitnessDistribution::Neutral => {
                let table = FitnessTable::from_model(sequence, n_symbols, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::Exponential(_) => {
                let table = FitnessTable::from_model(sequence, n_symbols, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::Lognormal(_) => {
                let table = FitnessTable::from_model(sequence, n_symbols, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::File(_) => {
                let table = FitnessTable::from_model(sequence, n_symbols, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::Epistatic(_, _) => {
                let table = FitnessTable::from_model(sequence, n_symbols, model)?;
                let epistasis = EpistasisTable::from_model(model)?;
                FitnessFunction::SimpleEpistatic(table, epistasis)
            }
        };
        Ok(Self {
            id,
            function,
            utility: model.utility.clone(),
        })
    }

    pub fn write(&self, path: &Path) -> Result<(), VirolutionError> {
        match &self.function {
            FitnessFunction::NonEpistatic(table) => {
                let table_name = format!("fitness_table_{}.npy", self.id);
                let mut table_file =
                    io::BufWriter::new(fs::File::create(path.join(table_name)).unwrap());
                table.write(&mut table_file).map_err(|e| {
                    VirolutionError::InitializationError(format!(
                        "Failed to write fitness table: {}",
                        e
                    ))
                })?;
            }
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                let table_name = format!("fitness_table_{}.npy", self.id);
                let mut table_file =
                    io::BufWriter::new(fs::File::create(path.join(table_name)).unwrap());
                table.write(&mut table_file).map_err(|e| {
                    VirolutionError::InitializationError(format!(
                        "Failed to write fitness table: {}",
                        e
                    ))
                })?;
                let epistasis_table_name = format!("fitness_table_{}.npy", self.id);
                let mut epistasis_table_file =
                    io::BufWriter::new(fs::File::create(path.join(epistasis_table_name)).unwrap());
                epistasis.write(&mut epistasis_table_file).map_err(|e| {
                    VirolutionError::InitializationError(format!(
                        "Failed to write epistasis table: {}",
                        e
                    ))
                })?;
            }
        }
        Ok(())
    }

    pub fn apply_utility(&self, fitness: f64) -> f64 {
        self.utility.apply(fitness)
    }

    pub fn get_fitness(&self, haplotype: &HaplotypeRef) -> f64 {
        match haplotype.try_get_changes() {
            Some(_) => {
                let fitness = haplotype
                    .try_get_ancestor()
                    .expect("Could not find ancestor during fitness update...")
                    .get_raw_fitness(self);
                let update = self.update_fitness(haplotype);
                fitness * update
            }
            None => self.compute_fitness(haplotype),
        }
    }

    fn update_fitness(&self, haplotype: &HaplotypeRef) -> f64 {
        match &self.function {
            FitnessFunction::NonEpistatic(table) => table.update_fitness(haplotype),
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                table.update_fitness(haplotype) * epistasis.update_fitness(haplotype)
            }
        }
    }

    fn compute_fitness(&self, haplotype: &HaplotypeRef) -> f64 {
        match &self.function {
            FitnessFunction::NonEpistatic(table) => table.compute_fitness(haplotype),
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                table.compute_fitness(haplotype) * epistasis.compute_fitness(haplotype)
            }
        }
    }
}
