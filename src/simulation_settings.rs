use super::fitness::FitnessDistribution;
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
    pub fitness_distribution: FitnessDistribution,
}

impl SimulationSettings {
    pub fn write(&self, writer: &mut dyn std::io::Write) {
        serde_yaml::to_writer(writer, self).expect("Unable to write settings.");
    }

    pub fn read(reader: &mut dyn std::io::Read) -> Self {
        serde_yaml::from_reader(reader).expect("Unable to read settings.")
    }

    pub fn write_to_file(&self, filename: &str) {
        let out =
            serde_yaml::to_string(self).unwrap_or_else(|_| panic!("Unable to serialize {self:?}."));
        fs::write(filename, out).unwrap_or_else(|_| panic!("Unable to write to {filename}."));
    }

    pub fn read_from_file(filename: &str) -> Self {
        let input = fs::read_to_string(filename)
            .unwrap_or_else(|_| panic!("Unable to read from {filename}."));
        serde_yaml::from_str(&input).unwrap_or_else(|_| panic!("Unable to deserialize {filename}."))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fitness::*;

    #[test]
    fn read_write_neutral() {
        let mut buffer = Vec::new();
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
            fitness_distribution: FitnessDistribution::Neutral,
        };
        settings.write(&mut buffer);
        let read_settings = SimulationSettings::read(&mut buffer.as_slice());
        assert_eq!(read_settings, settings);
    }

    #[test]
    fn read_write_exponential() {
        let mut buffer = Vec::new();
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
            fitness_distribution: FitnessDistribution::Exponential(ExponentialParameters {
                weights: MutationCategoryWeights {
                    beneficial: 0.29,
                    deleterious: 0.51,
                    lethal: 0.2,
                    neutral: 0.,
                },
                lambda_beneficial: 0.03,
                lambda_deleterious: 0.21,
            }),
        };
        settings.write(&mut buffer);
        let read_settings = SimulationSettings::read(&mut buffer.as_slice());
        assert_eq!(read_settings, settings);
    }

    #[test]
    fn read_write_file() {
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
            fitness_distribution: FitnessDistribution::Neutral,
        };
        settings.write_to_file("test_settings.yaml");
        let read_settings = SimulationSettings::read_from_file("test_settings.yaml");
        assert_eq!(read_settings, settings);
        std::fs::remove_file("test_settings.yaml").unwrap();
    }
}
