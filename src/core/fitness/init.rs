use npyz::NpyFile;
use rand::distributions::{Distribution, WeightedIndex};
use rand_distr::Exp;
use serde::{Deserialize, Serialize};

use crate::core::haplotype::Symbol;

use super::epistasis::EpiEntry;
use super::utility::UtilityFunction;

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
pub enum FitnessDistribution {
    Neutral,
    Exponential(ExponentialParameters),
    File(FileParameters),
    Lognormal(LognormalParameters),
    Epistatic(FileParameters, EpiFileParameters),
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

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct EpiFileParameters {
    pub path: String,
}

impl MutationCategoryWeights {
    fn as_slice(&self) -> &[f64] {
        unsafe { std::slice::from_raw_parts(self as *const Self as *const f64, 4) }
    }
}

impl ExponentialParameters {
    pub fn create_table(&self, n_symbols: usize, sequence: &[Symbol]) -> Vec<f64> {
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
    pub fn create_table(&self, n_symbols: usize, sequence: &[Symbol]) -> Vec<f64> {
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
    pub fn load_table(&self) -> Vec<f64> {
        let reader = NpyFile::new(std::fs::File::open(&self.path).unwrap()).unwrap();
        reader
            .data::<f64>()
            .unwrap()
            .map(|entry| entry.unwrap())
            .collect()
    }
}

impl EpiFileParameters {
    pub fn load_table(&self) -> Vec<EpiEntry> {
        let reader = NpyFile::new(std::fs::File::open(&self.path).unwrap()).unwrap();
        reader
            .data::<EpiEntry>()
            .unwrap()
            .map(|entry| entry.unwrap())
            .collect()
    }
}
