use std::fs;
use std::io;
use std::path::Path;

use crate::init::{FitnessDistribution, FitnessModel};

use crate::core::attributes::{AttributeProvider, AttributeProviderType, AttributeValue};
use crate::core::fitness::{EpistasisTable, FitnessTable, UtilityFunction};
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
    function: Box<dyn FitnessFunction<S>>,
    utility: UtilityFunction,
}

/// A trait for fitness functions.
pub trait FitnessFunction<S: Symbol>: Send + Sync + std::fmt::Debug {
    fn update_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64;
    fn compute_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64;
    fn write(&self, path: &Path, name: &'static str) -> Result<(), VirolutionError>;
    fn clone_box(&self) -> Box<dyn FitnessFunction<S>>;
}

impl<S: Symbol> Clone for Box<dyn FitnessFunction<S>> {
    fn clone(&self) -> Self {
        self.clone_box()
    }
}

#[derive(Clone, Debug)]
pub struct NonEpistatic<S: Symbol> {
    table: FitnessTable,
    phantom: std::marker::PhantomData<S>,
}

/// Non-epistatic fitness function implementation.
///
/// The fitness of a haplotype is computed from individual fitness values for each site in the
/// sequence. The computation is based on a precomputed fitness table.
impl<S: Symbol> NonEpistatic<S> {
    pub fn new(table: FitnessTable) -> Self {
        Self {
            table,
            phantom: std::marker::PhantomData,
        }
    }
}

impl<S: Symbol> FitnessFunction<S> for NonEpistatic<S> {
    fn update_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        self.table.update_fitness(haplotype)
    }

    fn compute_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        self.table.compute_fitness(haplotype)
    }

    fn write(&self, path: &Path, name: &'static str) -> Result<(), VirolutionError> {
        let table_name = format!("{name}_table.npy");
        let mut table_file = io::BufWriter::new(fs::File::create(path.join(table_name)).unwrap());
        self.table.write(&mut table_file)?;
        Ok(())
    }

    fn clone_box(&self) -> Box<dyn FitnessFunction<S>> {
        Box::new(self.clone())
    }
}

/// Simple epistatic fitness function implementation.
///
/// The fitness of a haplotype is computed from individual fitness values for each site in the
/// sequence. Additionally the fitness is modulated by an epistasis table which defines first order
/// interactions between sites.
#[derive(Clone, Debug)]
pub struct SimpleEpistatic<S: Symbol> {
    table: FitnessTable,
    epistasis: EpistasisTable<S>,
}

impl<S: Symbol> SimpleEpistatic<S> {
    pub fn new(table: FitnessTable, epistasis: EpistasisTable<S>) -> Self {
        Self { table, epistasis }
    }
}

impl<S: Symbol> FitnessFunction<S> for SimpleEpistatic<S> {
    fn update_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        self.table.update_fitness(haplotype) * self.epistasis.update_fitness(haplotype)
    }

    fn compute_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        self.table.compute_fitness(haplotype) * self.epistasis.compute_fitness(haplotype)
    }

    fn write(&self, path: &Path, name: &'static str) -> Result<(), VirolutionError> {
        let table_name = format!("{name}_table.npy");
        let epistasis_name = format!("{name}_epistasis_table.npy");

        let table_file = fs::File::create(path.join(table_name)).unwrap();
        let epistasis_path = fs::File::create(path.join(epistasis_name)).unwrap();

        let mut table_writer = io::BufWriter::new(table_file);
        let mut epistasis_writer = io::BufWriter::new(epistasis_path);

        self.table.write(&mut table_writer)?;
        self.epistasis.write(&mut epistasis_writer)?;
        Ok(())
    }

    fn clone_box(&self) -> Box<dyn FitnessFunction<S>> {
        Box::new(self.clone())
    }
}

/// Neutral fitness function implementation.
///
/// This avoids unnecessary lookups and computations for neutral fitness values.
#[derive(Clone, Debug)]
pub struct Neutral<S: Symbol> {
    marker: std::marker::PhantomData<S>,
}

impl<S: Symbol> Neutral<S> {
    pub fn new() -> Self {
        Self {
            marker: std::marker::PhantomData,
        }
    }
}

impl<S: Symbol> FitnessFunction<S> for Neutral<S> {
    fn update_fitness(&self, _haplotype: &HaplotypeRef<S>) -> f64 {
        1.0
    }

    fn compute_fitness(&self, _haplotype: &HaplotypeRef<S>) -> f64 {
        1.0
    }

    fn write(&self, _path: &Path, _name: &'static str) -> Result<(), VirolutionError> {
        Ok(())
    }

    fn clone_box(&self) -> Box<dyn FitnessFunction<S>> {
        Box::new(self.clone())
    }
}

impl<S: Symbol> FitnessProvider<S> {
    /// Create a new fitness provider.
    pub fn new(
        name: &'static str,
        function: Box<dyn FitnessFunction<S>>,
        utility: UtilityFunction,
    ) -> Self {
        Self {
            name,
            function,
            utility,
        }
    }

    /// Create a new fitness provider from a fitness model definition.
    ///
    /// Fitness models can be specified in the configuration file. Available fitness models are
    /// defined in `virolution::init`.
    pub fn from_model(
        name: &'static str,
        sequence: &[S],
        model: &FitnessModel,
    ) -> Result<Self, VirolutionError> {
        let function: Box<dyn FitnessFunction<S>> = match model.distribution {
            FitnessDistribution::Neutral => Box::new(Neutral::new()),
            FitnessDistribution::Exponential(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                Box::new(NonEpistatic::new(table))
            }
            FitnessDistribution::Lognormal(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                Box::new(NonEpistatic::new(table))
            }
            FitnessDistribution::File(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                Box::new(NonEpistatic::new(table))
            }
            FitnessDistribution::Epistatic(_) => {
                let table = FitnessTable::from_model(sequence, model)?;
                let epistasis = EpistasisTable::from_model(model)?;
                Box::new(SimpleEpistatic::new(table, epistasis))
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
                        .get_or_compute_attribute_raw(self.name)
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
        self.function.update_fitness(haplotype)
    }

    fn compute_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        self.function.compute_fitness(haplotype)
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
        self.function.write(path, self.name)
    }
}
