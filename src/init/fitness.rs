//! Fitness model and distribution definitions.
//!
//! Anything needed to initialize the [`crate::providers::FitnessProvider`] at the start of the
//! simulation.

use npyz::NpyFile;
use rand_distr::{weighted::WeightedIndex, Distribution, Exp};
use serde::{Deserialize, Serialize};

use crate::core::fitness::epistasis::EpiEntry;
use crate::core::fitness::UtilityFunction;
use crate::encoding::Symbol;

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

    /// Prepend a path to the file paths in the fitness model.
    ///
    /// This is used to make sure that paths can be defined relative to the configuration files
    /// location.
    pub(crate) fn prepend_path(&mut self, path: &str) {
        if path == "" {
            return;
        }
        match &mut self.distribution {
            FitnessDistribution::File(params) => {
                params.path = format!("{}/{}", path, params.path);
            }
            FitnessDistribution::Epistatic(params) => {
                params.path = format!("{}/{}", path, params.path);
                params.epi_path = format!("{}/{}", path, params.epi_path);
            }
            _ => {}
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum FitnessDistribution {
    Neutral,
    Exponential(ExponentialParameters),
    File(FileParameters),
    Lognormal(LognormalParameters),
    Epistatic(EpiFileParameters),
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
    pub epi_path: String,
}

impl MutationCategoryWeights {
    fn as_slice(&self) -> &[f64] {
        unsafe { std::slice::from_raw_parts(self as *const Self as *const f64, 4) }
    }
}

impl ExponentialParameters {
    pub fn create_table<S: Symbol>(&self, n_symbols: usize, sequence: &[S]) -> Vec<f64> {
        let mut table = vec![-1.; sequence.len() * n_symbols];
        let mut rng = rand::rng();

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
            table[position * n_symbols + symbol.index()] = 1.;
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
    pub fn create_table<S: Symbol>(&self, n_symbols: usize, sequence: &[S]) -> Vec<f64> {
        let mut table = vec![-1.; sequence.len() * n_symbols];
        let mut rng = rand::rng();

        // create distributions
        let lognormal = rand_distr::LogNormal::new(self.mu, self.sigma).unwrap();
        let lethal_distr = rand_distr::Bernoulli::new(self.lethal).unwrap();

        // fill wt sequence with ones
        for (position, &symbol) in sequence.iter().enumerate() {
            table[position * n_symbols + symbol.index()] = 1.;
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
        let reader = NpyFile::new(
            std::fs::File::open(&self.path)
                .unwrap_or_else(|_e| panic!("Could not open file: {}", self.path)),
        )
        .unwrap();
        reader
            .data::<f64>()
            .unwrap()
            .map(|entry| entry.unwrap())
            .collect()
    }
}

impl EpiFileParameters {
    pub fn load_table(&self) -> Vec<f64> {
        let file = std::fs::File::open(&self.path)
            .unwrap_or_else(|_e| panic!("Could not open file: {}", self.path));
        let reader = NpyFile::new(file).unwrap();
        reader
            .data::<f64>()
            .unwrap()
            .map(|entry| entry.unwrap())
            .collect()
    }

    pub fn load_epistasis(&self) -> Vec<EpiEntry> {
        let reader = NpyFile::new(std::fs::File::open(&self.epi_path).unwrap()).unwrap();
        reader
            .data::<EpiEntry>()
            .unwrap()
            .map(|entry| entry.unwrap())
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::fitness::FitnessTable;

    use crate::encoding::Nucleotide as Nt;

    #[test]
    fn create_neutral_table() {
        let nucleotides = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let sequence = vec![Nt::A; 100];
        let fitness = FitnessTable::from_model(
            &sequence,
            &FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
        )
        .unwrap();
        for position in 0..100 {
            for s in nucleotides.iter() {
                assert_eq!(fitness.get_value(&position, s), 1.);
            }
        }
    }

    #[test]
    fn create_exponential_table() {
        let nucleotides = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let sequence = vec![Nt::A; 100];
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
        for position in 0..100 {
            for s in nucleotides.iter() {
                if *s.index() == 0 {
                    assert_eq!(fitness.get_value(&position, s), 1.);
                } else {
                    assert_ne!(fitness.get_value(&position, s), 1.);
                }
            }
        }
    }
}
