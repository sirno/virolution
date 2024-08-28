use crate::core::haplotype::Symbol;
use crate::core::population::Population;
use crate::encoding::DECODE;
use crate::references::HaplotypeRef;
use serde::Deserialize;

pub trait PopulationIO {
    fn read(path: &str, wildtype: HaplotypeRef) -> Result<Population, csv::Error>;
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
    fn read(path: &str, wildtype: HaplotypeRef) -> Result<Population, csv::Error> {
        let mut reader = csv::Reader::from_path(path)?;
        let mut populations: Vec<Population> = Vec::new();

        for record in reader.deserialize() {
            let record: HaplotypeRecord = record?;
            if record.haplotype == "wt" {
                continue;
            }
            let mutations = record.haplotype.split(';');
            let mut positions: Vec<usize> = Vec::new();
            let mut changes: Vec<Symbol> = Vec::new();
            mutations.for_each(|mutation| {
                let mut mutation = mutation.split(':');
                let position = mutation.next().unwrap();

                let mut change = mutation.next().unwrap().split("->");
                let _origin = change.next();
                let target = change.next().unwrap().chars().next().unwrap() as u8;

                positions.push(position.parse::<usize>().unwrap());
                changes.push(DECODE.get(&target).copied());
            });
            let haplotype = wildtype.create_descendant(positions, changes, 0);
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
