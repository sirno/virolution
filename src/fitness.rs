use super::haplotype::Symbol;
use rand::prelude::*;
use rand_distr::{Exp, WeightedIndex};

#[derive(Clone)]
pub struct FitnessTable {
    // n_sites: usize,
    n_symbols: usize,
    table: Vec<f64>,
}

pub enum FitnessDistribution {
    Neutral,
    Exponential(ExponentialParameters),
    // Spikes,
    // Lognormal,
    // Beta,
    // Empirical,
}

pub enum MutationCategory {
    Beneficial,
    Deleterious,
    Lethal,
    Neutral,
}

pub struct MutationCategoryWeights {
    pub beneficial: f64,
    pub deleterious: f64,
    pub lethal: f64,
    pub neutral: f64,
}

pub struct ExponentialParameters {
    pub weights: MutationCategoryWeights,
    pub lambda_beneficial: f64,
    pub lambda_deleterious: f64,
}

impl MutationCategoryWeights {
    fn as_slice(&self) -> &[f64] {
        unsafe { std::slice::from_raw_parts(self as *const Self as *const f64, 4) }
    }
}

impl ExponentialParameters {
    fn create_table(&self, n_symbols: &usize, sequence: &Vec<Symbol>) -> Vec<f64> {
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
        for idx in 0..table.len() {
            if table[idx] >= 0. {
                continue;
            }

            match categories[category_distr.sample(&mut rng)] {
                MutationCategory::Beneficial => table[idx] = 1. + beneficial_distr.sample(&mut rng),
                MutationCategory::Deleterious => {
                    table[idx] = 1. - deleterious_distr.sample(&mut rng)
                }
                MutationCategory::Lethal => table[idx] = 0.,
                MutationCategory::Neutral => table[idx] = 1.,
            }
        }
        table
    }
}

impl FitnessTable {
    pub fn new(
        sequence: &Vec<Symbol>,
        n_symbols: &usize,
        distribution: FitnessDistribution,
    ) -> Self {
        let n_sites = sequence.len();
        let table = match distribution {
            FitnessDistribution::Exponential(params) => params.create_table(n_symbols, sequence),
            FitnessDistribution::Neutral => vec![1.; n_sites * n_symbols],
        };
        Self {
            // n_sites: n_sites,
            n_symbols: *n_symbols,
            table: table,
        }
    }

    pub fn get_fitness(&self, position: &usize, symbol: &Symbol) -> f64 {
        match symbol {
            Some(s) => self.table[position * self.n_symbols + *s as usize],
            None => 1.,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn create_neutral_table() {
        let sequence = vec![Some(0x00); 100];
        let fitness = FitnessTable::new(&sequence, &4, FitnessDistribution::Neutral);
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

        let fitness = FitnessTable::new(&sequence, &4, distribution);
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
        let fitness = FitnessTable::new(&sequence, &4, FitnessDistribution::Neutral);
        for position in 0..100 {
            for s in 0..4 {
                assert_eq!(fitness.get_fitness(&position, &Some(s)), 1.);
                assert_eq!(fitness.get_fitness(&position, &None), 1.);
            }
        }
    }
}
