use csv;
use evalexpr::context_map;
use phf::phf_map;
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;

use super::parameters::Parameters;

pub static TRANSFERS: phf::Map<&'static str, &[&[f64]]> = phf_map! {
    "migration_fwd" => MIGRATION_FWD,
    "migration_rev" => MIGRATION_REV,
    "root_ab" => ROOT_AB,
    "root_bc" => ROOT_BC,
    "root_cd" => ROOT_CD,
    "default" => DEFAULT,
};

const DEFAULT: &[&[f64]] = &[
    &[1., 0., 0., 0.],
    &[0., 1., 0., 0.],
    &[0., 0., 1., 0.],
    &[0., 0., 0., 1.],
];

const MIGRATION_FWD: &[&[f64]] = &[
    &[1., 0., 0., 0.],
    &[0.2, 0.8, 0., 0.],
    &[0., 0.2, 0.8, 0.],
    &[0., 0., 0., 0.],
];
const MIGRATION_REV: &[&[f64]] = &[
    &[0.8, 0.2, 0., 0.],
    &[0., 1., 0., 0.],
    &[0., 0., 1., 0.],
    &[0., 0., 0., 0.],
];

const ROOT_AB: &[&[f64]] = &[
    &[1., 0., 0., 0.],
    &[1., 1., 0., 0.],
    &[0., 0., 1., 0.],
    &[0., 0., 0., 1.],
];
const ROOT_BC: &[&[f64]] = &[
    &[1., 0., 0., 0.],
    &[0., 1., 0., 0.],
    &[0., 1., 1., 0.],
    &[0., 0., 0., 1.],
];
const ROOT_CD: &[&[f64]] = &[
    &[1., 0., 0., 0.],
    &[0., 1., 0., 0.],
    &[0., 0., 1., 0.],
    &[0., 0., 1., 1.],
];

#[derive(Clone, Debug, PartialEq)]
pub struct TransferMatrix<T> {
    matrix: Vec<T>,
    size: usize,
}

impl<T> TransferMatrix<T> {
    fn from_slice(slice: &[&[T]]) -> Self
    where
        T: Clone,
    {
        // Verify if slice is square
        slice
            .iter()
            .for_each(|row| assert_eq!(row.len(), slice.len()));

        let matrix = slice.iter().flat_map(|s| s.iter()).cloned().collect();
        Self {
            matrix,
            size: slice.len(),
        }
    }

    fn from_vec(vec: Vec<Vec<T>>) -> Self
    where
        T: Clone,
    {
        // Verify if vec is square
        vec.iter().for_each(|row| assert_eq!(row.len(), vec.len()));

        let matrix = vec.iter().flat_map(|s| s.iter()).cloned().collect();
        Self {
            matrix,
            size: vec.len(),
        }
    }

    pub fn get(&self, row: usize, col: usize) -> &T {
        &self.matrix[row * self.size + col]
    }

    pub fn get_row(&self, row: usize) -> &[T] {
        let begin = row * self.size;
        let end = begin + self.size;
        &self.matrix[begin..end]
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Schedule {
    table: Vec<ScheduleRecord>,
    transfers: HashMap<String, TransferMatrix<f64>>,
}

#[derive(Clone, Debug, PartialEq, Deserialize, Serialize)]
pub struct ScheduleRecord {
    generation: String,
    event: String,
    value: String,
}

#[derive(Debug)]
pub enum ScheduleError {
    IoError(std::io::Error),
    CsvError(csv::Error),
}

impl Serialize for Schedule {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(self.table.len()))?;
        for record in &self.table {
            seq.serialize_element(record)?;
        }
        seq.end()
    }
}

impl<'de> Deserialize<'de> for Schedule {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        let table: Vec<ScheduleRecord> = Vec::<ScheduleRecord>::deserialize(deserializer)?;
        Self::from_vec(table).map_err(|e| serde::de::Error::custom(format!("{:?}", e)))
    }
}

impl Schedule {
    pub fn read(filename: &str) -> Result<Self, ScheduleError> {
        let mut reader = BufReader::new(File::open(filename).map_err(ScheduleError::IoError)?);
        Schedule::from_reader(&mut reader)
    }

    pub fn from_reader(reader: &mut dyn std::io::Read) -> Result<Self, ScheduleError> {
        let mut reader = csv::ReaderBuilder::new()
            .delimiter(b';')
            .from_reader(reader);

        // Parse plan table
        let table: Vec<ScheduleRecord> = reader
            .deserialize()
            .collect::<Result<Vec<ScheduleRecord>, csv::Error>>()
            .map_err(ScheduleError::CsvError)?;

        Self::from_vec(table)
    }

    pub fn from_vec(table: Vec<ScheduleRecord>) -> Result<Self, ScheduleError> {
        let mut transfers: HashMap<String, TransferMatrix<f64>> = table
            .iter()
            .filter_map(|record| {
                if record.event != "transmission" {
                    return None;
                }

                if TRANSFERS.contains_key(record.value.as_str()) {
                    return Some((
                        record.value.clone(),
                        TransferMatrix::from_slice(TRANSFERS[record.value.as_str()]),
                    ));
                }

                let matrix: Vec<Vec<f64>> = serde_yaml::from_str(record.value.as_str()).unwrap();

                Some((record.value.clone(), TransferMatrix::from_vec(matrix)))
            })
            .collect();

        // Add default transfer matrix
        transfers.insert("default".to_string(), TransferMatrix::from_slice(DEFAULT));

        Ok(Self { table, transfers })
    }

    pub fn get_sample_size(&self, generation: usize) -> usize {
        match self.get_event_value("sample", generation) {
            Some(sample_size) => sample_size.parse().expect("Invalid value for sample"),
            None => 0,
        }
    }

    pub fn get_transfer_matrix(&self, generation: usize) -> &TransferMatrix<f64> {
        match self.get_event_value("transmission", generation) {
            Some(transfer_name) => &self.transfers[transfer_name],
            None => &self.transfers["default"],
        }
    }

    pub fn get_settings(&self, generation: usize) -> Option<Parameters> {
        self.get_event_value("settings", generation)
            .and_then(|settings_path| Parameters::read_from_file(settings_path).ok())
    }

    pub fn get_event_value(&self, event: &str, generation: usize) -> Option<&str> {
        match self
            .table
            .iter()
            .find(|record| match_generation(record, generation) && record.event == event)
        {
            Some(record) => Some(record.value.as_str()),
            None => None,
        }
    }

    pub fn check_transfer_table_sizes(&self, min: usize) -> bool {
        self.transfers.values().all(|t| t.size >= min)
    }
}

impl ScheduleRecord {
    pub fn new(generation: &str, event: &str, value: &str) -> Self {
        Self {
            generation: generation.to_string(),
            event: event.to_string(),
            value: value.to_string(),
        }
    }
}

fn match_generation(record: &ScheduleRecord, generation: usize) -> bool {
    match record.generation.parse::<usize>() {
        Ok(record_generation) => record_generation == generation,
        Err(_) => {
            let context = context_map! {
                "x" => generation as i64,
                "generation" => generation as i64,
                "{}" => generation as i64,
            }
            .unwrap();
            match evalexpr::eval_int_with_context(record.generation.as_str(), &context) {
                Ok(result) => result == 0,
                Err(_) => panic!("Invalid generation `{generation}`."),
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transfer_matrix() {
        let plan_content = r#"generation;event;value
1;transmission;migration_fwd
2;transmission;migration_rev
3;transmission;root_ab
4;transmission;root_bc
5;transmission;root_cd"#;

        let plan = Schedule::from_reader(&mut plan_content.as_bytes()).unwrap();

        assert_eq!(
            plan.get_transfer_matrix(0),
            &TransferMatrix::from_slice(TRANSFERS["default"])
        );
        assert_eq!(
            plan.get_transfer_matrix(1),
            &TransferMatrix::from_slice(TRANSFERS["migration_fwd"])
        );
        assert_eq!(
            plan.get_transfer_matrix(2),
            &TransferMatrix::from_slice(TRANSFERS["migration_rev"])
        );
        assert_eq!(
            plan.get_transfer_matrix(3),
            &TransferMatrix::from_slice(TRANSFERS["root_ab"])
        );
        assert_eq!(
            plan.get_transfer_matrix(4),
            &TransferMatrix::from_slice(TRANSFERS["root_bc"])
        );
        assert_eq!(
            plan.get_transfer_matrix(5),
            &TransferMatrix::from_slice(TRANSFERS["root_cd"])
        );
    }

    #[test]
    fn custom_transfer_matrix() {
        let plan_content = r#"generation;event;value
1;transmission;[[1, 2], [3, 4]]"#;

        let plan = Schedule::from_reader(&mut plan_content.as_bytes()).unwrap();
        assert_eq!(
            plan.get_transfer_matrix(1),
            &TransferMatrix::from_vec(vec![vec![1.0, 2.0], vec![3.0, 4.0]])
        );
    }

    #[test]
    fn sample_size() {
        let plan_content = r#"generation;event;value
1;sample;100
2;sample;200
3;sample;300"#;

        let plan = Schedule::from_reader(&mut plan_content.as_bytes()).unwrap();

        assert_eq!(plan.get_sample_size(0), 0);
        assert_eq!(plan.get_sample_size(1), 100);
        assert_eq!(plan.get_sample_size(2), 200);
        assert_eq!(plan.get_sample_size(3), 300);
    }

    #[test]
    fn expressions() {
        let plan_content = r#"generation;event;value
{} % 200;sample;100
{} % 10;transmission;migration_fwd
(5 + {}) % 10;transmission;migration_rev"#;

        let plan = Schedule::from_reader(&mut plan_content.as_bytes()).unwrap();

        for i in 0..=1000 {
            if i % 200 == 0 {
                assert_eq!(plan.get_sample_size(i), 100);
            } else {
                assert_eq!(plan.get_sample_size(i), 0);
            }

            if i % 10 == 0 {
                assert_eq!(
                    plan.get_transfer_matrix(i),
                    &TransferMatrix::from_slice(TRANSFERS["migration_fwd"])
                );
            } else if (5 + i) % 10 == 0 {
                assert_eq!(
                    plan.get_transfer_matrix(i),
                    &TransferMatrix::from_slice(TRANSFERS["migration_rev"])
                );
            } else {
                assert_eq!(
                    plan.get_transfer_matrix(i),
                    &TransferMatrix::from_slice(TRANSFERS["default"])
                );
            }
        }
    }
}
