use serde::{Deserialize, Serialize};
use std::fs;

use crate::init::FitnessModel;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Parameters {
    pub mutation_rate: f64,
    pub recombination_rate: f64,
    pub host_population_size: usize,
    pub infection_fraction: f64,
    pub basic_reproductive_number: f64,
    pub max_population: usize,
    pub dilution: f64,
    pub substitution_matrix: [[f64; 4]; 4],
    pub fitness_model: FitnessModelField,
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum FitnessModelField {
    SingleHost(FitnessModel),
    MultiHost(Vec<FitnessModelFraction>),
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct FitnessModelFraction {
    pub fraction: f64,
    pub fitness_model: FitnessModel,
}

#[derive(Debug)]
pub enum ParametersError {
    IoError(std::io::Error),
    YamlError(serde_yaml::Error),
}

impl std::fmt::Display for Parameters {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut output = vec![];
        self.write(&mut output).map_err(|_| std::fmt::Error)?;
        write!(formatter, "{}", String::from_utf8(output).unwrap())
    }
}

impl Parameters {
    pub fn write(&self, writer: &mut dyn std::io::Write) -> Result<(), ParametersError> {
        serde_yaml::to_writer(writer, self).map_err(ParametersError::YamlError)
    }

    pub fn read(reader: &mut dyn std::io::Read) -> Result<Parameters, ParametersError> {
        serde_yaml::from_reader(reader).map_err(ParametersError::YamlError)
    }

    pub fn write_to_file(&self, filename: &str) -> Result<(), ParametersError> {
        let file = fs::File::create(filename).map_err(ParametersError::IoError)?;
        let mut writer = std::io::BufWriter::new(file);
        self.write(&mut writer)
    }

    pub fn read_from_file(filename: &str) -> Result<Parameters, ParametersError> {
        let file = fs::File::open(filename).map_err(ParametersError::IoError)?;
        let mut reader = std::io::BufReader::new(file);
        Self::read(&mut reader)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::init::fitness::*;
    use crate::core::fitness::utility::*;

    #[test]
    fn read_write_neutral() {
        let mut buffer = Vec::new();
        let settings = Parameters {
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
            fitness_model: FitnessModelField::SingleHost(FitnessModel::new(
                FitnessDistribution::Neutral,
                UtilityFunction::Linear,
            )),
        };
        settings.write(&mut buffer).unwrap();
        let read_settings = Parameters::read(&mut buffer.as_slice()).unwrap();
        assert_eq!(read_settings, settings);
    }

    #[test]
    fn read_write_exponential() {
        let mut buffer = Vec::new();
        let settings = Parameters {
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
            fitness_model: FitnessModelField::SingleHost(FitnessModel::new(
                FitnessDistribution::Exponential(ExponentialParameters {
                    weights: MutationCategoryWeights {
                        beneficial: 0.29,
                        deleterious: 0.51,
                        lethal: 0.2,
                        neutral: 0.,
                    },
                    lambda_beneficial: 0.03,
                    lambda_deleterious: 0.21,
                }),
                UtilityFunction::Linear,
            )),
        };
        settings.write(&mut buffer).unwrap();
        let read_settings = Parameters::read(&mut buffer.as_slice()).unwrap();
        assert_eq!(read_settings, settings);
    }

    #[test]
    fn read_write_file() {
        let tmp_dir = std::env::temp_dir().join("test_settings.yaml");
        let path = tmp_dir.to_str().unwrap();
        let settings = Parameters {
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
            fitness_model: FitnessModelField::SingleHost(FitnessModel::new(
                FitnessDistribution::Neutral,
                UtilityFunction::Linear,
            )),
        };
        settings.write_to_file(path).unwrap();
        let read_settings = Parameters::read_from_file(path).unwrap();
        assert_eq!(read_settings, settings);
        std::fs::remove_file(path).unwrap();
    }
}
