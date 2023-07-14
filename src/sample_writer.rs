use csv;
use seq_io::fasta::{OwnedRecord, Record};
use std::{collections::HashMap, fs, io, path::Path};

use itertools::Itertools;

use crate::{barcode::BarcodeEntry, haplotype::FASTA_ENCODE, simulation::SimulationTrait};

pub trait SampleWriter {
    fn write(&self, simulations: &[Box<SimulationTrait>], sample_size: usize);
}

pub struct FastaSampleWriter<'a> {
    simulation_name: &'a str,
    path: &'a str,
}

impl<'a> FastaSampleWriter<'a> {
    pub fn new(simulation_name: &'a str, path: &'a str) -> Self {
        // init barcode file header
        let barcode_path = Path::new(path).join("barcodes.csv");
        let mut barcode_file = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .open(barcode_path)
            .expect("Unable to open barcode file.");
        BarcodeEntry::write_header(&mut barcode_file)
            .expect("Unable to write header to barcode file.");

        Self {
            simulation_name,
            path,
        }
    }
}

pub struct CsvSampleWriter<'a> {
    simulation_name: &'a str,
    path: &'a str,
}

impl<'a> CsvSampleWriter<'a> {
    pub fn new(simulation_name: &'a str, path: &'a str) -> Self {
        // init barcode file header
        let barcode_path = Path::new(path).join("barcodes.csv");
        let mut barcode_file = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .open(barcode_path)
            .expect("Unable to open barcode file.");
        BarcodeEntry::write_header(&mut barcode_file)
            .expect("Unable to write header to barcode file.");

        Self {
            simulation_name,
            path,
        }
    }
}

impl<'a> SampleWriter for FastaSampleWriter<'a> {
    fn write(&self, simulations: &[Box<SimulationTrait>], sample_size: usize) {
        log::info!("Writing fasta sample to {}", self.path);
        for (compartment_id, compartment) in simulations.iter().enumerate() {
            let generation = compartment.get_generation();
            let barcode = format!("sample_{generation}_{compartment_id}");

            // create output files
            let barcode_path = Path::new(self.path).join("barcodes.csv");
            let sample_path = Path::new(self.path).join(format!("{barcode}.fasta"));

            // create output directories
            std::fs::create_dir_all(barcode_path.parent().unwrap()).unwrap_or_else(|_| {
                eprintln!("Unable to create output path.");
                std::process::exit(1);
            });

            // create file buffers
            let mut barcode_file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(barcode_path)
                .expect("Unable to open barcode file.");
            let mut samples_file = io::BufWriter::new(fs::File::create(sample_path).unwrap());

            // precompute sequences
            let population = compartment.get_population();
            let sequences: HashMap<usize, Vec<u8>> = population
                .iter()
                .unique()
                .map(|haplotype_ref| {
                    let sequence = haplotype_ref
                        .get_sequence()
                        .into_iter()
                        .map(|symbol| match symbol {
                            Some(s) => FASTA_ENCODE[&s],
                            None => 0x2d,
                        })
                        .collect();
                    (haplotype_ref.get_id(), sequence)
                })
                .collect();

            // sample sequences and write to file
            for (haplotype_id, haplotype_ref) in population
                .choose_multiple(&mut rand::thread_rng(), sample_size)
                .into_iter()
                .enumerate()
            {
                let head = format!(
                    "compartment_id={};sequence_id={};generation={}",
                    compartment_id, haplotype_id, generation
                )
                .as_bytes()
                .to_vec();
                let sequence = sequences[&haplotype_ref.get_id()].clone();
                let record = OwnedRecord {
                    head,
                    seq: sequence,
                };
                record
                    .write(&mut samples_file)
                    .expect("Unable to write to file.");
            }

            // write barcode to file
            BarcodeEntry {
                barcode: &barcode,
                experiment: &self.simulation_name.to_string(),
                time: generation,
                replicate: 0,
                compartment: compartment_id,
            }
            .write(&mut barcode_file)
            .expect("Unable to write to barcode file.");
        }
    }
}

impl<'a> SampleWriter for CsvSampleWriter<'a> {
    fn write(&self, simulations: &[Box<SimulationTrait>], sample_size: usize) {
        log::info!("Writing csv sample to {}", self.path);
        for (compartment_id, compartment) in simulations.iter().enumerate() {
            let generation = compartment.get_generation();
            let barcode = format!("sample_{generation}_{compartment_id}");

            // create output files
            let barcode_path = Path::new(self.path).join("barcodes.csv");
            let sample_path = Path::new(self.path).join(format!("{barcode}.csv"));

            // create output directories
            std::fs::create_dir_all(barcode_path.parent().unwrap()).unwrap_or_else(|_| {
                eprintln!("Unable to create output path.");
                std::process::exit(1);
            });

            // create file buffers
            let mut barcode_file = fs::OpenOptions::new()
                .create(true)
                .append(true)
                .open(barcode_path)
                .expect("Unable to open barcode file.");
            let mut samples_file = csv::WriterBuilder::new()
                .from_path(sample_path)
                .expect("Unable to open sample file.");

            samples_file
                .write_record(&["haplotype", "count"])
                .expect("Unable to write header to samples file.");

            // sample sequences and write to file
            let population = compartment.get_population();
            for (haplotype_ref, haplotype_count) in population
                .choose_multiple(&mut rand::thread_rng(), sample_size)
                .into_iter()
                .counts()
            {
                samples_file
                    .write_record(&[
                        haplotype_ref.as_ref().get_string(),
                        haplotype_count.to_string(),
                    ])
                    .expect("Unable to write to samples file.")
            }

            // write barcode to file
            BarcodeEntry {
                barcode: &barcode,
                experiment: &self.simulation_name.to_string(),
                time: generation,
                replicate: 0,
                compartment: compartment_id,
            }
            .write(&mut barcode_file)
            .expect("Unable to write to barcode file.");
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_csv_write() {
        let tmp_dir = std::env::temp_dir();
        let path = tmp_dir.to_str().unwrap();
        let _sample_writer = CsvSampleWriter::new("test_simulation", path);
    }
}
