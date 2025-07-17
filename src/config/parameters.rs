use serde::{Deserialize, Serialize};
use std::cmp::min;
use std::fs;
use std::path::Path;
use std::sync::Arc;

use crate::core::AttributeSetDefinition;
use crate::core::hosts::{HostSpec, HostSpecs};
use crate::encoding::Symbol;
use crate::init::FitnessModel;
use crate::providers::FitnessProvider;
use crate::simulation::BasicHost;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Parameters {
    /// The mutation rate represents the probability of a mutation occurring in a virus at each
    /// time step.
    pub mutation_rate: f64,

    /// The recombination rate represents the probability of recombination occurring between two
    /// viruses in the population.
    pub recombination_rate: f64,

    /// The size of the host population in the simulation.
    pub host_population_size: usize,

    /// The basic reproductive number represents the average burst size of a virus in the population.
    pub basic_reproductive_number: f64,

    /// The maximum population size for number of viruses in the simulation.
    pub max_population: usize,

    /// The dilution factor represents the fraction of the population that is diluted at each
    /// time step.
    pub dilution: f64,

    /// The substitution matrix is a 4x4 matrix representing the probabilities of nucleotide
    /// substitutions.
    pub substitution_matrix: [[f64; 4]; 4],

    /// The host model defines host population structure and fitness characteristics.
    pub host_model: HostModel,

    /// The number of times a virus can hit hosts before it fails to infect.
    #[serde(default = "default_n_hits")]
    pub n_hits: usize,
}

fn default_n_hits() -> usize {
    1
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct HostFitness {
    pub reproductive: Option<FitnessModel>,
    pub infective: Option<FitnessModel>,
}

impl HostFitness {
    pub fn new(reproductive: Option<FitnessModel>, infective: Option<FitnessModel>) -> Self {
        Self {
            reproductive,
            infective,
        }
    }

    pub fn make_definitions<S: Symbol>(
        &self,
        parameters: &Parameters,
        sequence: &[S],
        path: Option<&Path>,
        host_id: usize,
    ) -> (AttributeSetDefinition<S>, BasicHost) {
        let mut attributes = AttributeSetDefinition::new();
        let path = path.unwrap_or_else(|| Path::new("./"));
        let reproductive_fitness_name = self.reproductive.as_ref().map(|model| {
            let mut model = model.clone();
            let name: &'static str =
                Box::leak(format!("reproductive_fitness_{host_id}").into_boxed_str());
            model.prepend_path(path.to_str().unwrap());
            attributes.register(Arc::new(
                FitnessProvider::from_model(name, sequence, &model).unwrap(),
            ));
            name
        });
        let infective_fitness_name = self.infective.as_ref().map(|model| {
            let mut model = model.clone();
            let name: &'static str =
                Box::leak(format!("infective_fitness_{host_id}").into_boxed_str());
            model.prepend_path(path.to_str().unwrap());
            attributes.register(Arc::new(
                FitnessProvider::from_model(name, sequence, &model).unwrap(),
            ));
            name
        });
        let host = BasicHost::new(
            sequence.len(),
            parameters,
            reproductive_fitness_name,
            infective_fitness_name,
        );
        (attributes, host)
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub enum HostModel {
    SingleHost(HostFitness),
    MultiHost(Vec<FitnessModelFraction>),
}

impl HostModel {
    pub fn make_definitions<S: Symbol>(
        &self,
        parameters: &Parameters,
        sequence: &[S],
        path: Option<&Path>,
    ) -> (AttributeSetDefinition<S>, HostSpecs<S>) {
        match self {
            HostModel::SingleHost(host_fitness) => {
                let (attributes, host) =
                    host_fitness.make_definitions(parameters, sequence, path, 0);
                let host_spec = HostSpec::new(0..parameters.host_population_size, Box::new(host));
                (attributes, HostSpecs::from_vec(vec![host_spec]))
            }
            HostModel::MultiHost(models) => {
                let mut attributes = AttributeSetDefinition::new();
                let mut hosts = HostSpecs::new();
                for (i, host_fraction) in models.iter().enumerate() {
                    let (attrs, host) = host_fraction
                        .host_fitness
                        .make_definitions(parameters, sequence, path, i);
                    let n_hosts = (host_fraction.fraction * parameters.host_population_size as f64)
                        .round() as usize;
                    let lower = i * n_hosts;
                    let upper = min(lower + n_hosts, parameters.host_population_size);
                    let host_spec = HostSpec::new(lower..upper, Box::new(host));
                    attributes.extend(attrs.into_iter());
                    hosts.push(host_spec);
                }
                (attributes, hosts)
            }
        }
    }
}

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct FitnessModelFraction {
    pub fraction: f64,
    pub host_fitness: HostFitness,
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
    use crate::core::fitness::utility::*;
    use crate::init::fitness::*;

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
            basic_reproductive_number: 100.,
            max_population: 100,
            dilution: 0.17,
            host_model: HostModel::SingleHost(HostFitness::new(
                Some(FitnessModel::new(
                    FitnessDistribution::Neutral,
                    UtilityFunction::Linear,
                )),
                None,
            )),
            n_hits: 1,
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
            basic_reproductive_number: 100.,
            max_population: 100,
            dilution: 0.17,
            host_model: HostModel::SingleHost(HostFitness::new(
                Some(FitnessModel::new(
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
                None,
            )),
            n_hits: 1,
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
            basic_reproductive_number: 100.,
            max_population: 100,
            dilution: 0.17,
            host_model: HostModel::SingleHost(HostFitness::new(
                Some(FitnessModel::new(
                    FitnessDistribution::Neutral,
                    UtilityFunction::Linear,
                )),
                None,
            )),
            n_hits: 1,
        };
        settings.write_to_file(path).unwrap();
        let read_settings = Parameters::read_from_file(path).unwrap();
        assert_eq!(read_settings, settings);
        std::fs::remove_file(path).unwrap();
    }
}
