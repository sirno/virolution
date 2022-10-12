use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SimulationSettings {
    pub mutation_rate: f64,
    pub recombination_rate: f64,
    pub host_population_size: usize,
    pub infection_fraction: f64,
    pub basic_reproductive_number: f64,
    pub max_population: usize,
    pub dilution: f64,
    pub substitution_matrix: [[f64; 4]; 4],
}

impl SimulationSettings {
    pub fn write(&self, filename: &str) {
        let out = serde_yaml::to_string(self)
            .unwrap_or_else(|_| panic!("Unable to serialize {:?}.", self));
        fs::write(filename, out).unwrap_or_else(|_| panic!("Unable to write to {}.", filename));
    }

    pub fn read(filename: &str) -> Self {
        let input = fs::read_to_string(filename)
            .unwrap_or_else(|_| panic!("Unable to read from {}.", filename));
        serde_yaml::from_str(&input)
            .unwrap_or_else(|_| panic!("Unable to deserialize {}.", filename))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn read_write() {
        let settings = SimulationSettings {
            mutation_rate: 1e-6,
            recombination_rate: 1e-8,
            substitution_matrix: [
                [0., 1., 1., 1.],
                [1., 0., 1., 1.],
                [1., 1., 0., 1.],
                [1., 1., 1., 0.],
            ],
            host_population_size: 5,
            infection_fraction: 0.7,
            basic_reproductive_number: 100.,
            max_population: 100,
            dilution: 0.17,
        };
        settings.write("test.yaml");
        let read_settings = SimulationSettings::read("test.yaml");
        assert_eq!(read_settings, settings);
        fs::remove_file("test.yaml").expect("Unable to remove file.");
    }
}
