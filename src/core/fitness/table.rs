use npyz::WriterBuilder;

use crate::init::{FitnessDistribution, FitnessModel};

use crate::encoding::Symbol;
use crate::errors::VirolutionError;
use crate::references::HaplotypeRef;

#[derive(Clone, Debug)]
pub struct FitnessTable {
    n_sites: usize,
    n_symbols: usize,
    table: Vec<f64>,
}

impl FitnessTable {
    pub fn new(n_sites: usize, n_symbols: usize, table: Vec<f64>) -> Self {
        Self {
            n_sites,
            n_symbols,
            table,
        }
    }

    pub fn from_model<S: Symbol>(
        sequence: &[S],
        fitness_model: &FitnessModel,
    ) -> Result<Self, VirolutionError> {
        let n_sites = sequence.len();
        let table = match fitness_model.distribution {
            FitnessDistribution::Neutral => vec![1.; n_sites * S::SIZE],
            FitnessDistribution::Exponential(ref params) => params.create_table(S::SIZE, sequence),
            FitnessDistribution::Lognormal(ref params) => params.create_table(S::SIZE, sequence),
            FitnessDistribution::File(ref params) => params.load_table(),
            FitnessDistribution::Epistatic(ref params) => params.load_table(),
        };

        if table.len() != n_sites * S::SIZE {
            return Err(VirolutionError::InitializationError(format!(
                "Fitness table has wrong size: {} instead of {}",
                table.len(),
                n_sites * S::SIZE
            )));
        }

        Ok(Self::new(n_sites, S::SIZE, table))
    }

    /// Compute factor to update the fitness of a haplotype change based on the nonepistatic table
    pub fn update_fitness<S: Symbol>(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        // Callers should ensure that the haplotype is mutant
        if !haplotype.is_mutant() {
            panic!("Cannot update fitness of haplotype that is not mutant");
        }

        // If there are not changes then there is no update
        let changes = match haplotype.try_get_changes() {
            Some(changes) => changes,
            None => {
                return 1.;
            }
        };

        // Compute the fitness update
        changes.iter().fold(1., |acc, change| {
            acc * self.get_value(&change.position, &change.to)
                / self.get_value(&change.position, &change.from)
        })
    }

    /// Compute the fitness of a haplotype without any epistatic effects
    pub fn compute_fitness<S: Symbol>(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        haplotype
            .get_mutations()
            .iter()
            .fold(1., |acc, (pos, symbol)| {
                acc * self.get_value(pos, &S::decode(symbol))
            })
    }

    #[inline]
    pub fn get_value<S: Symbol>(&self, position: &usize, symbol: &S) -> f64 {
        self.table[*position * S::SIZE + *symbol.index()]
    }

    pub fn write(&self, writer: &mut impl std::io::Write) -> Result<(), VirolutionError> {
        let shape = &[self.n_sites as u64, self.n_symbols as u64];
        let mut npy_writer = npyz::WriteOptions::new()
            .default_dtype()
            .shape(shape)
            .writer(writer)
            .begin_nd()
            .map_err(|e| VirolutionError::WriteError(format!("{}", e)))?;
        npy_writer
            .extend(self.table.clone())
            .map_err(|e| VirolutionError::WriteError(format!("{}", e)))?;
        npy_writer
            .finish()
            .map_err(|e| VirolutionError::WriteError(format!("{}", e)))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::UtilityFunction;
    use super::*;
    use crate::encoding::Nucleotide as Nt;
    use crate::init::fitness::{ExponentialParameters, MutationCategoryWeights};

    #[test]
    fn get_value() {
        let symbols = &[Nt::A, Nt::T, Nt::C, Nt::G];
        let sequence = vec![Nt::A; 100];
        let fitness = FitnessTable::from_model(
            &sequence,
            &FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
        )
        .unwrap();
        for position in 0..100 {
            for s in symbols.iter() {
                assert_eq!(fitness.get_value(&position, s), 1.);
            }
        }
    }

    #[test]
    fn get_value_indexing() {
        let table = FitnessTable {
            n_sites: 2,
            n_symbols: 4,
            table: vec![1., 2., 3., 4., 5., 6., 7., 8.],
        };
        assert_eq!(table.get_value(&0, &Nt::A), 1.);
        assert_eq!(table.get_value(&0, &Nt::T), 2.);
        assert_eq!(table.get_value(&0, &Nt::C), 3.);
        assert_eq!(table.get_value(&0, &Nt::G), 4.);
        assert_eq!(table.get_value(&1, &Nt::A), 5.);
        assert_eq!(table.get_value(&1, &Nt::T), 6.);
        assert_eq!(table.get_value(&1, &Nt::C), 7.);
        assert_eq!(table.get_value(&1, &Nt::G), 8.);
    }

    #[test]
    fn write_table() {
        let sequence = vec![Nt::A; 100000];
        let distribution = FitnessDistribution::Exponential(ExponentialParameters {
            weights: MutationCategoryWeights {
                beneficial: 0.29,
                deleterious: 0.51,
                lethal: 0.2,
                neutral: 0.,
            },
            lambda_beneficial: 0.03,
            lambda_deleterious: 0.21,
        });

        let fitness = FitnessTable::from_model(
            &sequence,
            &FitnessModel::new(distribution, UtilityFunction::Linear),
        )
        .unwrap();
        let mut buffer = Vec::new();

        fitness.write(&mut buffer).unwrap();

        let npy_data = npyz::NpyFile::new(buffer.as_slice()).unwrap();
        assert_eq!(npy_data.shape(), &[100000, 4]);
        let data: Vec<f64> = npy_data
            .data::<f64>()
            .unwrap()
            .map(|el| el.unwrap())
            .collect();

        assert_eq!(data, fitness.table);
    }
}
