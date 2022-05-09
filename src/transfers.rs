use csv;
use phf::phf_map;
use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;

pub static TRANSFERS: phf::Map<&'static str, &[[f64; 4]; 4]> = phf_map! {
    "migration_fwd" => &MIGRATION_FWD,
    "migration_rev" => &MIGRATION_REV,
};

const MIGRATION: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [0., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 0.],
];

const MIGRATION_FWD: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [0.2, 0.8, 0., 0.],
    [0., 0.2, 0.8, 0.],
    [0., 0., 0., 0.],
];
const MIGRATION_REV: [[f64; 4]; 4] = [
    [0.8, 0.2, 0., 0.],
    [0., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 0.],
];

#[derive(Debug)]
pub struct Plan {
    table: Vec<PlanRecord>,
}

#[derive(Debug, Deserialize)]
struct PlanRecord {
    generation: usize,
    event: String,
    value: String,
}

impl Plan {
    pub fn read(filename: &str) -> Self {
        let file = File::open(filename).expect("File not found.");
        let mut reader = csv::Reader::from_reader(BufReader::new(file));
        let table: Vec<PlanRecord> = reader
            .deserialize()
            .map(|record| record.expect("Unable to deserialize record."))
            .collect();
        Self { table: table }
    }

    pub fn get_transfer_matrix(&self, generation: usize) -> &[[f64; 4]; 4] {
        match self.get_transfer_name(generation) {
            Some(transfer_name) => TRANSFERS[transfer_name],
            None => &MIGRATION,
        }
    }

    #[inline]
    pub fn get_transfer_name(&self, generation: usize) -> Option<&str> {
        match self
            .table
            .iter()
            .find(|record| record.generation == generation && record.event == "transmission")
        {
            Some(record) => Some(record.value.as_str()),
            None => None,
        }
    }
}
