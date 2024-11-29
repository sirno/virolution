use std::fs;
use std::io;
use std::path::Path;

use super::epistasis::EpistasisTable;
use super::init::{FitnessDistribution, FitnessModel};
use super::table::FitnessTable;
use super::utility::UtilityFunction;

use crate::core::attributes::{AttributeProvider, AttributeProviderType, AttributeValue};
use crate::encoding::Symbol;
use crate::errors::VirolutionError;
use crate::references::HaplotypeRef;

/// A provider of fitness values for haplotypes.
///
/// The fitness provider is responsible for computing the fitness of a haplotype based on its
/// sequence and the fitness model.
#[derive(Clone, Debug)]
pub struct FitnessProvider<S: Symbol> {
    name: &'static str,
    function: FitnessFunction<S>,
    utility: UtilityFunction,
}

#[derive(Clone, Debug)]
pub enum FitnessFunction<S: Symbol> {
    NonEpistatic(FitnessTable),
    SimpleEpistatic(FitnessTable, EpistasisTable<S>),
}

impl<S: Symbol> FitnessProvider<S> {
    /// Create a new fitness provider from a fitness model.
    ///
    /// Fitness models can be specified in the configuration file. Available fitness models are
    /// defined in `virolution::core::fitness::init`.
    pub fn from_model(
        name: &'static str,
        sequence: &[S],
        model: &FitnessModel,
    ) -> Result<Self, VirolutionError> {
        let function = match model.distribution {
            FitnessDistribution::Neutral => {
                let table = FitnessTable::from_model(sequence, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::Exponential(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::Lognormal(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::File(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                FitnessFunction::NonEpistatic(table)
            }
            FitnessDistribution::Epistatic(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                let epistasis = EpistasisTable::from_model(model)?;
                FitnessFunction::SimpleEpistatic(table, epistasis)
            }
        };
        Ok(Self {
            name,
            function,
            utility: model.utility.clone(),
        })
    }

    pub fn apply_utility(&self, fitness: f64) -> f64 {
        self.utility.apply(fitness)
    }

    pub fn get_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        match haplotype.try_get_changes() {
            Some(_) => {
                let fitness = f64::try_from(
                    haplotype
                        .try_get_ancestor()
                        .expect("Could not find ancestor during fitness update...")
                        .get_attribute_or_compute(self.name)
                        .unwrap(),
                )
                .unwrap();
                let update = self.update_fitness(haplotype);
                fitness * update
            }
            None => self.compute_fitness(haplotype),
        }
    }

    fn update_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        match &self.function {
            FitnessFunction::NonEpistatic(table) => table.update_fitness(haplotype),
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                table.update_fitness(haplotype) * epistasis.update_fitness(haplotype)
            }
        }
    }

    fn compute_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        match &self.function {
            FitnessFunction::NonEpistatic(table) => table.compute_fitness(haplotype),
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                table.compute_fitness(haplotype) * epistasis.compute_fitness(haplotype)
            }
        }
    }
}

impl<S: Symbol> AttributeProvider<S> for FitnessProvider<S> {
    fn name(&self) -> &'static str {
        self.name
    }

    fn get_provider_type(&self) -> AttributeProviderType {
        AttributeProviderType::Lazy
    }

    fn compute(&self, haplotype: &Option<HaplotypeRef<S>>) -> AttributeValue {
        match haplotype {
            Some(haplotype) => AttributeValue::F64(self.get_fitness(haplotype)),
            None => panic!("Haplotype not found during fitness computation."),
        }
    }

    fn map(&self, value: AttributeValue) -> AttributeValue {
        match value {
            AttributeValue::F64(fitness) => AttributeValue::F64(self.apply_utility(fitness)),
            _ => value,
        }
    }

    fn write(&self, path: &Path) -> Result<(), VirolutionError> {
        match &self.function {
            FitnessFunction::NonEpistatic(table) => {
                let table_name = format!("{}_table.npy", self.name);
                let mut table_file =
                    io::BufWriter::new(fs::File::create(path.join(table_name)).unwrap());
                table.write(&mut table_file)?;
            }
            FitnessFunction::SimpleEpistatic(table, epistasis) => {
                let table_name = format!("{}_table.npy", self.name);
                let epistasis_name = format!("{}_epistasis_table.npy", self.name);

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
}
