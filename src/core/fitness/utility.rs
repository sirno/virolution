use serde::{Deserialize, Deserializer, Serialize};

/// Utility function used to transform fitness values into parameter values.
#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum UtilityFunction {
    /// Linear utility function.
    Linear,

    /// Algebraic utility function with an upper bound `upper` to avoid infinite values and/or get
    /// diminishing returns.
    Algebraic { upper: f64 },

    /// Hill utility function with parameters `p` and `n`. The function maps fitness values to
    /// values between 0 and 1. Values from (0, 1) are mapped to (0, p) and values from (1, ∞) are
    /// mapped to (p, 1).
    #[serde(deserialize_with = "deserialize_hill")]
    Hill {
        p: f64,
        n: f64,
        #[serde(skip_serializing)]
        k: f64,
    },
}

// Custom deserialization function for the Hill variant
fn deserialize_hill<'de, D>(deserializer: D) -> Result<(f64, f64, f64), D::Error>
where
    D: Deserializer<'de>,
{
    // Use a temporary struct to deserialize p and n
    #[derive(Deserialize)]
    struct TempHill {
        p: f64,
        n: f64,
    }

    let temp = TempHill::deserialize(deserializer)?;

    if temp.n < 0.0 || temp.p > 1.0 {
        return Err(serde::de::Error::custom("p must be between 0 and 1"));
    }

    let k = (1.0 - temp.p) / temp.p;
    Ok((temp.p, temp.n, k))
}

impl UtilityFunction {
    pub fn apply(&self, fitness: f64) -> f64 {
        match self {
            UtilityFunction::Linear => fitness,
            UtilityFunction::Algebraic { upper } => {
                let factor = 1. / (upper - 1.);
                upper * factor * fitness / (1. + factor * fitness)
            }
            UtilityFunction::Hill { n, k, .. } => {
                let fitness_pow_n = fitness.powf(*n);
                fitness_pow_n / (k + fitness_pow_n)
            }
        }
    }
}
