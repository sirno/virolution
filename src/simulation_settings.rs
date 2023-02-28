use super::simulation_parameters::SimulationParameters;
use super::simulation_plan::SimulationPlan;
use serde::{Deserialize, Serialize};
use std::fs;

#[derive(Serialize, Deserialize, Debug, Clone, PartialEq)]
pub struct SimulationSettings {
    pub simulation_settings: SimulationParameters,
    pub simulation_plan: SimulationPlan,
}

#[derive(Debug)]
pub enum SimulationSettingsError {
    IoError(std::io::Error),
    YamlError(serde_yaml::Error),
}

impl SimulationSettings {
    pub fn write(&self, writer: &mut dyn std::io::Write) -> Result<(), SimulationSettingsError> {
        serde_yaml::to_writer(writer, self).map_err(SimulationSettingsError::YamlError)
    }

    pub fn read(
        reader: &mut dyn std::io::Read,
    ) -> Result<SimulationSettings, SimulationSettingsError> {
        serde_yaml::from_reader(reader).map_err(SimulationSettingsError::YamlError)
    }

    pub fn write_to_file(&self, filename: &str) -> Result<(), SimulationSettingsError> {
        let file = fs::File::create(filename).map_err(SimulationSettingsError::IoError)?;
        let mut writer = std::io::BufWriter::new(file);
        self.write(&mut writer)
    }

    pub fn read_from_file(filename: &str) -> Result<SimulationSettings, SimulationSettingsError> {
        let file = fs::File::open(filename).map_err(SimulationSettingsError::IoError)?;
        let mut reader = std::io::BufReader::new(file);
        Self::read(&mut reader)
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::simulation_plan::SimulationPlanRecord;

    #[test]
    fn read_write() {
        let settings = SimulationSettings {
            simulation_settings: SimulationParameters {
                mutation_rate: 0.1,
                recombination_rate: 0.1,
                host_population_size: 100,
                infection_fraction: 0.1,
                basic_reproductive_number: 1.0,
                max_population: 1000,
                dilution: 0.1,
                substitution_matrix: [[0.0; 4]; 4],
                fitness_model: crate::fitness::FitnessModel {
                    distribution: crate::fitness::FitnessDistribution::Neutral,
                    utility: crate::fitness::UtilityFunction::Linear,
                },
            },
            simulation_plan: SimulationPlan::from_vec(vec![
                SimulationPlanRecord::new("0", "transmission", "migration_fwd"),
                SimulationPlanRecord::new("1", "transmission", "migration_rev"),
                SimulationPlanRecord::new("2", "sample", "1000"),
            ])
            .unwrap(),
        };
        let mut output = vec![];
        settings.write(&mut output).unwrap();
        let settings2 = SimulationSettings::read(&mut &output[..]).unwrap();
        assert_eq!(settings, settings2);
    }
}
