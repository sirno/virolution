use super::haplotype::Symbol;
use npyz::{NpyFile, WriterBuilder};
use rand::prelude::*;
use rand_distr::{Exp, WeightedIndex};
use serde::{Deserialize, Serialize};

#[derive(Clone)]
pub struct FitnessTable {
    n_sites: usize,
    n_symbols: usize,
    table: Vec<f64>,
    fitness_model: FitnessModel,
}

#[derive(Clone, Debug)]
pub struct FitnessError {
    pub message: String,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum FitnessDistribution {
    Neutral,
    Exponential(ExponentialParameters),
    File(FileParameters),
    Lognormal(LognormalParameters),
    // Spikes,
    // Beta,
    // Empirical,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum MutationCategory {
    Beneficial,
    Deleterious,
    Lethal,
    Neutral,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct MutationCategoryWeights {
    pub beneficial: f64,
    pub deleterious: f64,
    pub lethal: f64,
    pub neutral: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct ExponentialParameters {
    pub weights: MutationCategoryWeights,
    pub lambda_beneficial: f64,
    pub lambda_deleterious: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct LognormalParameters {
    pub lethal: f64,
    pub mu: f64,
    pub sigma: f64,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct FileParameters {
    pub path: String,
}

impl MutationCategoryWeights {
    fn as_slice(&self) -> &[f64] {
        unsafe { std::slice::from_raw_parts(self as *const Self as *const f64, 4) }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct FitnessModel {
    pub distribution: FitnessDistribution,
    pub utility: UtilityFunction,
}

impl FitnessModel {
    pub fn new(distribution: FitnessDistribution, utility: UtilityFunction) -> Self {
        Self {
            distribution,
            utility,
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum UtilityFunction {
    Linear,
    Algebraic { upper: f64 },
}

impl ExponentialParameters {
    fn create_table(&self, n_symbols: usize, sequence: &Vec<Symbol>) -> Vec<f64> {
        let mut table = vec![-1.; sequence.len() * n_symbols];
        let mut rng = rand::thread_rng();

        // create distributions
        let categories = [
            MutationCategory::Beneficial,
            MutationCategory::Deleterious,
            MutationCategory::Lethal,
            MutationCategory::Neutral,
        ];
        let weights = self.weights.as_slice();
        let category_distr = WeightedIndex::new(weights).unwrap();
        let beneficial_distr = Exp::new(1. / self.lambda_beneficial).unwrap();
        let deleterious_distr = Exp::new(1. / self.lambda_deleterious).unwrap();

        // fill wt sequence with ones
        for (position, &symbol) in sequence.iter().enumerate() {
            if let Some(value) = symbol {
                table[position * n_symbols + (value as usize)] = 1.;
            }
        }

        // sample for remaining positions
        for el in &mut table {
            if *el >= 0. {
                continue;
            }

            match categories[category_distr.sample(&mut rng)] {
                MutationCategory::Beneficial => *el = 1. + beneficial_distr.sample(&mut rng),
                MutationCategory::Deleterious => {
                    *el = f64::max(1. - deleterious_distr.sample(&mut rng), 0.)
                }
                MutationCategory::Lethal => *el = 0.,
                MutationCategory::Neutral => *el = 1.,
            }
        }

        table
    }
}

impl LognormalParameters {
    fn create_table(&self, n_symbols: usize, sequence: &Vec<Symbol>) -> Vec<f64> {
        let mut table = vec![-1.; sequence.len() * n_symbols];
        let mut rng = rand::thread_rng();

        // create distributions
        let lognormal = rand_distr::LogNormal::new(self.mu, self.sigma).unwrap();
        let lethal_distr = rand_distr::Bernoulli::new(self.lethal).unwrap();

        // fill wt sequence with ones
        for (position, &symbol) in sequence.iter().enumerate() {
            if let Some(value) = symbol {
                table[position * n_symbols + (value as usize)] = 1.;
            }
        }

        // sample for remaining positions
        for el in &mut table {
            if *el >= 0. {
                continue;
            }

            *el = if lethal_distr.sample(&mut rng) {
                0.
            } else {
                lognormal.sample(&mut rng)
            };
        }

        table
    }
}

impl FileParameters {
    fn load_table(&self) -> Vec<f64> {
        let reader = NpyFile::new(std::fs::File::open(&self.path).unwrap()).unwrap();
        reader
            .data::<f64>()
            .unwrap()
            .map(|entry| entry.unwrap())
            .collect()
    }
}

impl FitnessTable {
    pub fn new(
        n_sites: usize,
        n_symbols: usize,
        table: Vec<f64>,
        fitness_model: FitnessModel,
    ) -> Self {
        Self {
            n_sites,
            n_symbols,
            table,
            fitness_model,
        }
    }

    pub fn from_model(
        sequence: &Vec<Symbol>,
        n_symbols: usize,
        fitness_model: FitnessModel,
    ) -> Result<Self, FitnessError> {
        let n_sites = sequence.len();
        let table = match fitness_model.distribution {
            FitnessDistribution::Exponential(ref params) => {
                params.create_table(n_symbols, sequence)
            }
            FitnessDistribution::Lognormal(ref params) => params.create_table(n_symbols, sequence),
            FitnessDistribution::File(ref params) => params.load_table(),
            FitnessDistribution::Neutral => vec![1.; n_sites * n_symbols],
        };

        if table.len() != n_sites * n_symbols {
            return Err(FitnessError {
                message: format!(
                    "Fitness table has wrong size: {} instead of {}",
                    table.len(),
                    n_sites * n_symbols
                ),
            });
        }

        Ok(Self::new(n_sites, n_symbols, table, fitness_model))
    }

    pub fn get_fitness(&self, position: &usize, symbol: &Symbol) -> f64 {
        match symbol {
            Some(s) => self.table[position * self.n_symbols + *s as usize],
            None => 1.,
        }
    }

    pub fn write(&self, writer: &mut impl std::io::Write) -> Result<(), std::io::Error> {
        let shape = &[self.n_sites as u64, self.n_symbols as u64];
        let mut npy_writer = npyz::WriteOptions::new()
            .default_dtype()
            .shape(shape)
            .writer(writer)
            .begin_nd()?;
        npy_writer.extend(self.table.clone())?;
        npy_writer.finish()?;
        Ok(())
    }

    pub fn utility(&self, fitness: f64) -> f64 {
        match self.fitness_model.utility {
            UtilityFunction::Linear => fitness,
            UtilityFunction::Algebraic { upper } => {
                let factor = 1. / (upper - 1.);
                upper * factor * fitness / (1. + factor * fitness)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_neutral_table() {
        let sequence = vec![Some(0x00); 100];
        let fitness = FitnessTable::from_model(
            &sequence,
            4,
            FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
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
            FitnessModel::new(distribution, UtilityFunction::Linear),
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
    fn get_fitness() {
        let sequence = vec![Some(0x00); 100];
        let fitness = FitnessTable::from_model(
            &sequence,
            4,
            FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
        )
        .unwrap();
        for position in 0..100 {
            for s in 0..4 {
                assert_eq!(fitness.get_fitness(&position, &Some(s)), 1.);
                assert_eq!(fitness.get_fitness(&position, &None), 1.);
            }
        }
    }

    #[test]
    fn get_fitness_indexing() {
        let table = FitnessTable {
            n_sites: 2,
            n_symbols: 4,
            table: vec![1., 2., 3., 4., 5., 6., 7., 8.],
            fitness_model: FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
        };
        assert_eq!(table.get_fitness(&0, &Some(0)), 1.);
        assert_eq!(table.get_fitness(&0, &Some(1)), 2.);
        assert_eq!(table.get_fitness(&0, &Some(2)), 3.);
        assert_eq!(table.get_fitness(&0, &Some(3)), 4.);
        assert_eq!(table.get_fitness(&1, &Some(0)), 5.);
        assert_eq!(table.get_fitness(&1, &Some(1)), 6.);
        assert_eq!(table.get_fitness(&1, &Some(2)), 7.);
        assert_eq!(table.get_fitness(&1, &Some(3)), 8.);
    }

    #[test]
    fn write_fitness() {
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
            FitnessModel::new(distribution, UtilityFunction::Linear),
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
