use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Serialize, Deserialize, Debug)]
pub struct SimulationSettings {
    pub mutation_rate: f64,
    pub substitution_matrix: [[f64; 4]; 4],
    pub host_population_size: usize,
    pub infection_fraction: f64,
    pub basic_reproductive_number: f64,
    pub max_population: usize,
    pub dilution: f64,
}

impl SimulationSettings {
    pub fn write(&self, filename: &str) {
        let out =
            serde_yaml::to_string(self).expect(format!("Unable to serialize {:?}.", self).as_str());
        fs::write(filename, out).expect(format!("Unable to write to {}.", filename).as_str());
    }

    pub fn read(filename: &str) -> Self {
        let input = fs::read_to_string(filename)
            .expect(format!("Unable to read from {}.", filename).as_str());
        serde_yaml::from_str(&input).expect(format!("Unable to deserialize {}.", filename).as_str())
    }
}
