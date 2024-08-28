use npyz::WriterBuilder;

use super::init::{FitnessDistribution, FitnessModel};

use crate::core::haplotype::Symbol;
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

    pub fn from_model(
        sequence: &[Symbol],
        n_symbols: usize,
        fitness_model: &FitnessModel,
    ) -> Result<Self, VirolutionError> {
        let n_sites = sequence.len();
        let table = match fitness_model.distribution {
            FitnessDistribution::Exponential(ref params) => {
                params.create_table(n_symbols, sequence)
            }
            FitnessDistribution::Lognormal(ref params) => params.create_table(n_symbols, sequence),
            FitnessDistribution::File(ref params) => params.load_table(),
            FitnessDistribution::Neutral => vec![1.; n_sites * n_symbols],
            FitnessDistribution::Epistatic(_, _) => todo!(),
        };

        if table.len() != n_sites * n_symbols {
            return Err(VirolutionError::InitializationError(format!(
                "Fitness table has wrong size: {} instead of {}",
                table.len(),
                n_sites * n_symbols
            )));
        }

        Ok(Self::new(n_sites, n_symbols, table))
    }

    /// Compute factor to update the fitness of a haplotype change based on the nonepistatic table
    pub fn update_fitness(&self, haplotype: &HaplotypeRef) -> f64 {
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
        changes.iter().fold(1., |acc, (position, (old, new))| {
            acc * self.get_value(&position, new) / self.get_value(&position, old)
        })
    }

    /// Compute the fitness of a haplotype without any epistatic effects
    pub fn compute_fitness(&self, haplotype: &HaplotypeRef) -> f64 {
        haplotype
            .get_mutations()
            .iter()
            .fold(1., |acc, (pos, symbol)| acc * self.get_value(&pos, symbol))
    }

    #[inline]
    pub fn get_value(&self, position: &usize, symbol: &Symbol) -> f64 {
        match symbol {
            Some(s) => self.table[position * self.n_symbols + *s as usize],
            None => 1.,
        }
    }

    pub fn write(&self, writer: &mut impl std::io::Write) -> Result<(), VirolutionError> {
        let shape = &[self.n_sites as u64, self.n_symbols as u64];
        let mut npy_writer = npyz::WriteOptions::new()
            .default_dtype()
            .shape(shape)
            .writer(writer)
            .begin_nd()
            .map_err(|e| VirolutionError::InitializationError(format!("{}", e)))?;
        npy_writer
            .extend(self.table.clone())
            .map_err(|e| VirolutionError::InitializationError(format!("{}", e)))?;
        npy_writer
            .finish()
            .map_err(|e| VirolutionError::InitializationError(format!("{}", e)))?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::super::init::{ExponentialParameters, MutationCategoryWeights};
    use super::super::utility::UtilityFunction;
    use super::*;

    #[test]
    fn create_neutral_table() {
        let sequence = vec![Some(0x00); 100];
        let fitness = FitnessTable::from_model(
            &sequence,
            4,
            &FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
        )
        .unwrap();
        for val in fitness.table {
            assert_eq!(val, 1.);
        }
    }

    #[test]
    fn create_exponential_table() {
        let sequence = vec![Some(0x00); 100];
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
            4,
            &FitnessModel::new(distribution, UtilityFunction::Linear),
        )
        .unwrap();
        for (idx, val) in fitness.table.into_iter().enumerate() {
            if idx % 4 == 0 {
                assert_eq!(val, 1.);
            } else {
                assert_ne!(val, 1.);
            }
        }
    }
    #[test]
    fn get_value() {
        let sequence = vec![Some(0x00); 100];
        let fitness = FitnessTable::from_model(
            &sequence,
            4,
            &FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
        )
        .unwrap();
        for position in 0..100 {
            for s in 0..4 {
                assert_eq!(fitness.get_value(&position, &Some(s)), 1.);
                assert_eq!(fitness.get_value(&position, &None), 1.);
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
        assert_eq!(table.get_value(&0, &Some(0)), 1.);
        assert_eq!(table.get_value(&0, &Some(1)), 2.);
        assert_eq!(table.get_value(&0, &Some(2)), 3.);
        assert_eq!(table.get_value(&0, &Some(3)), 4.);
        assert_eq!(table.get_value(&1, &Some(0)), 5.);
        assert_eq!(table.get_value(&1, &Some(1)), 6.);
        assert_eq!(table.get_value(&1, &Some(2)), 7.);
        assert_eq!(table.get_value(&1, &Some(3)), 8.);
    }

    #[test]
    fn write_table() {
        let sequence = vec![Some(0x00); 100000];
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
            4,
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
