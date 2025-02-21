//! Settings module.

use super::parameters::Parameters;
use super::schedule::Schedule;

use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct Settings {
    pub parameters: Vec<Parameters>,
    pub schedule: Schedule,
}

#[derive(Debug)]
pub enum SettingsError {
    IoError(std::io::Error),
    YamlError(serde_yaml::Error),
}

impl std::error::Error for SettingsError {}

impl std::fmt::Display for SettingsError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SettingsError::IoError(error) => write!(formatter, "IO error: {}", error),
            SettingsError::YamlError(error) => write!(formatter, "YAML error: {}", error),
        }
    }
}

impl std::fmt::Display for Settings {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut output = vec![];
        self.write(&mut output).map_err(|_| std::fmt::Error)?;
        write!(formatter, "{}", String::from_utf8(output).unwrap())
    }
}

impl Settings {
    pub fn write(&self, writer: &mut dyn std::io::Write) -> Result<(), SettingsError> {
        serde_yaml::to_writer(writer, self).map_err(SettingsError::YamlError)
    }

    pub fn read(reader: &mut dyn std::io::Read) -> Result<Settings, SettingsError> {
        serde_yaml::from_reader(reader).map_err(SettingsError::YamlError)
    }

    pub fn write_to_file(&self, filename: &str) -> Result<(), SettingsError> {
        let file = fs::File::create(filename).map_err(SettingsError::IoError)?;
        let mut writer = std::io::BufWriter::new(file);
        self.write(&mut writer)
    }

    pub fn read_from_file(filename: &str) -> Result<Settings, SettingsError> {
        let file = fs::File::open(filename).map_err(SettingsError::IoError)?;
        let mut reader = std::io::BufReader::new(file);
        Self::read(&mut reader)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::config::parameters::FitnessModelField;
    use crate::config::schedule::ScheduleRecord;
    use crate::init::FitnessModel;

    #[test]
    fn read_write() {
        let settings = Settings {
            parameters: vec![Parameters {
                mutation_rate: 0.1,
                recombination_rate: 0.1,
                host_population_size: 100,
                infection_fraction: 0.1,
                basic_reproductive_number: 1.0,
                max_population: 1000,
                dilution: 0.1,
                substitution_matrix: [[0.0; 4]; 4],
                fitness_model: FitnessModelField::SingleHost(FitnessModel::new(
                    crate::init::fitness::FitnessDistribution::Neutral,
                    crate::core::fitness::utility::UtilityFunction::Linear,
                )),
            }],
            schedule: Schedule::from_vec(vec![
                ScheduleRecord::new("0", "transmission", "migration_fwd"),
                ScheduleRecord::new("1", "transmission", "migration_rev"),
                ScheduleRecord::new("2", "sample", "1000"),
            ])
            .unwrap(),
        };
        let mut output = vec![];
        settings.write(&mut output).unwrap();
        let settings2 = Settings::read(&mut &output[..]).unwrap();
        assert_eq!(settings, settings2);
    }
}
