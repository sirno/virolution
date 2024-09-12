use serde::Deserialize;

use crate::core::haplotype::Symbol;
use crate::core::population::Population;
use crate::encoding::DECODE;
use crate::errors::{Result, VirolutionError};
use crate::references::HaplotypeRef;

pub trait PopulationIO {
    fn read(path: &str, wildtype: HaplotypeRef) -> Result<Population>;
    fn write(&self, path: &str);
}

#[derive(Debug, Deserialize)]
struct HaplotypeRecord {
    haplotype: String,
    count: usize,
}

impl PopulationIO for Population {
    /// Reads a CSV file containing haplotypes and their counts.
    ///
    /// Warning: This function will create multiple instances of the same
    /// haplotype if it is present multiple times in the CSV file.
    fn read(path: &str, wildtype: HaplotypeRef) -> Result<Population> {
        let mut reader = csv::Reader::from_path(path)
            .map_err(|_err| VirolutionError::ReadError(format!("Failed to read from {path}")))?;
        let mut populations: Vec<Population> = Vec::new();

        for record in reader.deserialize() {
            let record: HaplotypeRecord = record.map_err(|_err| {
                VirolutionError::ReadError(format!("Failed to parse record in {path}"))
            })?;

            // skip wildtype
            if record.haplotype == "wt" {
                continue;
            }

            // parse haplotype
            let (positions, changes) = parse_haplotype(&record.haplotype)?;
            let haplotype = wildtype.create_descendant(positions, changes, 0);

            // create population and add for merging
            let population = Population::from_haplotype(haplotype, record.count);
            populations.push(population);
        }
        let population = Population::from_iter(populations);
        Ok(population)
    }

    fn write(&self, _path: &str) {
        unimplemented!()
    }
}

fn parse_haplotype(haplotype: &str) -> Result<(Vec<usize>, Vec<Symbol>)> {
    let mutations = haplotype.split(';');
    let parsed_mutations: Result<Vec<(usize, Symbol)>> = mutations.map(parse_mutation).collect();
    Ok(parsed_mutations?.into_iter().unzip())
}

fn parse_mutation(mutation: &str) -> Result<(usize, Symbol)> {
    let mutation_split: Vec<&str> = mutation.split(':').collect();
    match mutation_split.as_slice() {
        [position, change] => Ok((
            position.parse::<usize>().map_err(|_| {
                VirolutionError::ReadError(format!("Failed to parse position: {position}"))
            })?,
            parse_change(change)?,
        )),
        _ => Err(VirolutionError::ReadError(format!(
            "Invalid mutation format: {mutation}"
        ))),
    }
}

fn parse_change(change: &str) -> Result<Symbol> {
    let change_split: Vec<&str> = change.split(':').collect();
    match change_split.as_slice() {
        [_origin, target] => {
            let target = target.chars().next().unwrap() as u8;
            Ok(DECODE.get(&target).copied())
        }
        _ => Err(VirolutionError::ReadError(format!(
            "Invalid change format: {change}"
        ))),
    }
}
