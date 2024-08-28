use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum UtilityFunction {
    Linear,
    Algebraic { upper: f64 },
}

impl UtilityFunction {
    pub fn apply(&self, fitness: f64) -> f64 {
        match self {
            UtilityFunction::Linear => fitness,
            UtilityFunction::Algebraic { upper } => {
                let factor = 1. / (upper - 1.);
                upper * factor * fitness / (1. + factor * fitness)
            }
        }
    }
}
