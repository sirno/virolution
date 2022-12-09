use super::fitness::FitnessModel;
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
    pub fitness_model: FitnessModel,
}

#[derive(Debug)]
pub enum SimulationSettingsError {
    IoError(std::io::Error),
    YamlError(serde_yaml::Error),
}

impl std::fmt::Display for SimulationSettings {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut output = vec![];
        self.write(&mut output).map_err(|_| std::fmt::Error)?;
        write!(formatter, "{}", String::from_utf8(output).unwrap())
    }
}

impl SimulationSettings {
    pub fn write(&self, writer: &mut dyn std::io::Write) -> Result<(), SimulationSettingsError> {
        serde_yaml::to_writer(writer, self).map_err(|err| SimulationSettingsError::YamlError(err))
    }

    pub fn read(
        reader: &mut dyn std::io::Read,
    ) -> Result<SimulationSettings, SimulationSettingsError> {
        serde_yaml::from_reader(reader).map_err(|err| SimulationSettingsError::YamlError(err))
    }

    pub fn write_to_file(&self, filename: &str) -> Result<(), SimulationSettingsError> {
        let file =
            fs::File::create(filename).map_err(|err| SimulationSettingsError::IoError(err))?;
        let mut writer = std::io::BufWriter::new(file);
        self.write(&mut writer)
    }

    pub fn read_from_file(filename: &str) -> Result<SimulationSettings, SimulationSettingsError> {
        let file = fs::File::open(filename).map_err(|err| SimulationSettingsError::IoError(err))?;
        let mut reader = std::io::BufReader::new(file);
        Self::read(&mut reader)
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
            fitness_model: FitnessModel {
                distribution: FitnessDistribution::Neutral,
                utility: UtilityFunction::Linear,
            },
        };
        settings.write(&mut buffer).unwrap();
        let read_settings = SimulationSettings::read(&mut buffer.as_slice()).unwrap();
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
            fitness_model: FitnessModel {
                distribution: FitnessDistribution::Exponential(ExponentialParameters {
                    weights: MutationCategoryWeights {
                        beneficial: 0.29,
                        deleterious: 0.51,
                        lethal: 0.2,
                        neutral: 0.,
                    },
                    lambda_beneficial: 0.03,
                    lambda_deleterious: 0.21,
                }),
                utility: UtilityFunction::Linear,
            },
        };
        settings.write(&mut buffer).unwrap();
        let read_settings = SimulationSettings::read(&mut buffer.as_slice()).unwrap();
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
            fitness_model: FitnessModel {
                distribution: FitnessDistribution::Neutral,
                utility: UtilityFunction::Linear,
            },
        };
        settings.write_to_file("test_settings.yaml").unwrap();
        let read_settings = SimulationSettings::read_from_file("test_settings.yaml").unwrap();
        assert_eq!(read_settings, settings);
        std::fs::remove_file("test_settings.yaml").unwrap();
    }
}
