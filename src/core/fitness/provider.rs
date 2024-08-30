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

/// A provider of fitness values for haplotypes.
///
/// The fitness provider is responsible for computing the fitness of a haplotype based on its
/// sequence and the fitness model.
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
    /// Create a new fitness provider from a fitness model.
    ///
    /// Fitness models can be specified in the configuration file. Available fitness models are
    /// defined in `virolution::core::fitness::init`.
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
            FitnessDistribution::Epistatic(_) => {
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
                table.write(&mut table_file)?;
            }
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                let table_name = format!("fitness_table_{}.npy", self.id);
                let epistasis_name = format!("fitness_table_{}.npy", self.id);

                let table_file = fs::File::create(path.join(table_name)).unwrap();
                let epistasis_path = fs::File::create(path.join(epistasis_name)).unwrap();

                let mut table_writer = io::BufWriter::new(table_file);
                let mut epistasis_writer = io::BufWriter::new(epistasis_path);

                table.write(&mut table_writer)?;
                epistasis.write(&mut epistasis_writer)?;
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
