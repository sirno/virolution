use csv;
use derive_more::{Deref, DerefMut};
use evalexpr::context_map;
use phf::phf_map;
use serde::Deserialize;
use std::fs::File;
use std::io::BufReader;

pub static TRANSFERS: phf::Map<&'static str, &[[f64; 4]; 4]> = phf_map! {
    "migration_fwd" => &MIGRATION_FWD,
    "migration_rev" => &MIGRATION_REV,
    "root_ab" => &ROOT_AB,
    "root_bc" => &ROOT_BC,
    "root_cd" => &ROOT_CD,
};

const DEFAULT: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [0., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 1.],
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

const ROOT_AB: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [1., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 0., 1.],
];
const ROOT_BC: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [0., 1., 0., 0.],
    [0., 1., 1., 0.],
    [0., 0., 0., 1.],
];
const ROOT_CD: [[f64; 4]; 4] = [
    [1., 0., 0., 0.],
    [0., 1., 0., 0.],
    [0., 0., 1., 0.],
    [0., 0., 1., 1.],
];

#[derive(Debug, Deref, DerefMut)]
pub struct Plan(Vec<PlanRecord>);

#[derive(Debug, Deserialize)]
pub struct PlanRecord {
    generation: String,
    event: String,
    value: String,
}

impl Plan {
    pub fn read(filename: &str) -> Self {
        let mut reader = BufReader::new(File::open(filename).expect("File not found."));
        Plan::from_reader(&mut reader)
    }

    pub fn from_reader(reader: &mut dyn std::io::Read) -> Self {
        let mut reader = csv::Reader::from_reader(reader);
        let table: Vec<PlanRecord> = reader
            .deserialize()
            .map(|record| record.expect("Unable to deserialize record."))
            .collect();
        Self(table)
    }

    pub fn get_sample_size(&self, generation: usize) -> usize {
        match self
            .iter()
            .find(|record| match_generation(record, generation) && record.event == "sample")
        {
            Some(record) => record.value.parse().expect("Invalid value for sample"),
            None => 0,
        }
    }

    pub fn get_transfer_matrix(&self, generation: usize) -> &[[f64; 4]; 4] {
        match self.get_transfer_name(generation) {
            Some(transfer_name) => TRANSFERS[transfer_name],
            None => &DEFAULT,
        }
    }

    #[inline]
    pub fn get_transfer_name(&self, generation: usize) -> Option<&str> {
        match self
            .iter()
            .find(|record| match_generation(record, generation) && record.event == "transmission")
        {
            Some(record) => Some(record.value.as_str()),
            None => None,
        }
    }
}

fn match_generation(record: &PlanRecord, generation: usize) -> bool {
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
        let plan_content = r#"generation,event,value
1,transmission,migration_fwd
2,transmission,migration_rev
3,transmission,root_ab
4,transmission,root_bc
5,transmission,root_cd"#;

        let plan = Plan::from_reader(&mut plan_content.as_bytes());

        assert_eq!(plan.get_transfer_matrix(0), &DEFAULT);
        assert_eq!(plan.get_transfer_matrix(1), &MIGRATION_FWD);
        assert_eq!(plan.get_transfer_matrix(2), &MIGRATION_REV);
        assert_eq!(plan.get_transfer_matrix(3), &ROOT_AB);
        assert_eq!(plan.get_transfer_matrix(4), &ROOT_BC);
        assert_eq!(plan.get_transfer_matrix(5), &ROOT_CD);
    }

    #[test]
    fn sample_size() {
        let plan_content = r#"generation,event,value
1,sample,100
2,sample,200
3,sample,300"#;

        let plan = Plan::from_reader(&mut plan_content.as_bytes());

        assert_eq!(plan.get_sample_size(0), 0);
        assert_eq!(plan.get_sample_size(1), 100);
        assert_eq!(plan.get_sample_size(2), 200);
        assert_eq!(plan.get_sample_size(3), 300);
    }

    #[test]
    fn expressions() {
        let plan_content = r#"generation,event,value
{} % 200,sample,100
{} % 10,transmission,migration_fwd
(5 + {}) % 10,transmission,migration_rev"#;

        let plan = Plan::from_reader(&mut plan_content.as_bytes());

        for i in 0..=1000 {
            if i % 200 == 0 {
                assert_eq!(plan.get_sample_size(i), 100);
            } else {
                assert_eq!(plan.get_sample_size(i), 0);
            }

            if i % 10 == 0 {
                assert_eq!(plan.get_transfer_matrix(i), &MIGRATION_FWD);
            } else if (5 + i) % 10 == 0 {
                assert_eq!(plan.get_transfer_matrix(i), &MIGRATION_REV);
            } else {
                assert_eq!(plan.get_transfer_matrix(i), &DEFAULT);
            }
        }
    }
}
