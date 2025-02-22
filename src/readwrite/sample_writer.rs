use seq_io::fasta::{OwnedRecord, Record};
use std::{
    cell::RefCell,
    collections::HashMap,
    fs, io,
    path::{Path, PathBuf},
    rc::Rc,
};

use itertools::Itertools;

use crate::core::population::Store;
use crate::core::{Historian, Population};
use crate::encoding::Symbol;
use crate::references::HaplotypeRefTrait;
use crate::{barcode::BarcodeEntry, simulation::SimulationTrait};

pub trait SampleWriter<S: Symbol + 'static> {
    /// Get the name of the experiment
    fn get_experiment_name(&self) -> &str;

    /// Get the path to the output directory
    fn get_path(&self) -> PathBuf;

    /// Get the historian
    fn get_historian(&self) -> Option<Rc<RefCell<Historian>>>;

    /// Write a sample of each simulation to a file
    fn sample_from_simulations(
        &self,
        simulations: &[Box<SimulationTrait<S>>],
        sample_size: usize,
    ) -> Result<(), std::io::Error> {
        log::info!("Writing sample to {}", self.get_path().to_str().unwrap());
        for (compartment_id, compartment) in simulations.iter().enumerate() {
            let generation = compartment.get_generation();

            // sample population
            let population = compartment.get_population();
            let sample = population.choose_multiple(&mut rand::rng(), sample_size);

            // write to file
            let barcode = self.write(&sample, compartment.get_generation(), compartment_id)?;

            if let Some(historian) = &self.get_historian() {
                historian
                    .borrow_mut()
                    .record_sample(generation, compartment_id, sample.clone());
            }

            // write barcode to file
            self.write_barcode(&BarcodeEntry {
                barcode: &barcode,
                experiment: &self.get_experiment_name().to_string(),
                time: generation,
                replicate: 0,
                compartment: compartment_id,
            })?;
        }
        Ok(())
    }

    /// Write a sample to a file and return the barcode
    fn write(
        &self,
        population: &Population<Store<S>>,
        generation: usize,
        compartment: usize,
    ) -> Result<String, std::io::Error>;

    /// Get the path to the barcode file
    fn get_barcode_path(&self) -> std::path::PathBuf;

    /// Create the barcode file and write the header
    fn init_barcode_file(&self) -> Result<(), std::io::Error> {
        let barcode_path = self.get_barcode_path();
        std::fs::create_dir_all(barcode_path.parent().unwrap())?;
        let mut barcode_file = fs::OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(barcode_path)?;
        BarcodeEntry::write_header(&mut barcode_file)?;
        Ok(())
    }

    /// Append a barcode to the barcode file
    fn write_barcode(&self, barcode: &BarcodeEntry) -> Result<(), std::io::Error> {
        let barcode_path = self.get_barcode_path();
        let mut barcode_file = fs::OpenOptions::new().append(true).open(barcode_path)?;
        barcode.write(&mut barcode_file)?;
        Ok(())
    }
}

pub struct FastaSampleWriter<'a> {
    simulation_name: &'a str,
    path: &'a str,
    historian: Option<Rc<RefCell<Historian>>>,
}

impl<'a> FastaSampleWriter<'a> {
    pub fn new<S: Symbol + 'static>(
        simulation_name: &'a str,
        path: &'a str,
        historian: Option<Rc<RefCell<Historian>>>,
    ) -> Result<Self, std::io::Error> {
        let writer = Self {
            simulation_name,
            path,
            historian,
        };

        SampleWriter::<S>::init_barcode_file(&writer)?;

        Ok(writer)
    }
}

pub struct CsvSampleWriter<'a> {
    simulation_name: &'a str,
    path: &'a str,
    historian: Option<Rc<RefCell<Historian>>>,
}

impl<'a> CsvSampleWriter<'a> {
    pub fn new<S: Symbol + 'static>(
        simulation_name: &'a str,
        path: &'a str,
        historian: Option<Rc<RefCell<Historian>>>,
    ) -> Result<Self, std::io::Error> {
        let writer = Self {
            simulation_name,
            path,
            historian,
        };

        SampleWriter::<S>::init_barcode_file(&writer)?;

        Ok(writer)
    }
}

impl<S: Symbol + 'static> SampleWriter<S> for FastaSampleWriter<'_> {
    fn get_experiment_name(&self) -> &str {
        self.simulation_name
    }

    fn get_path(&self) -> PathBuf {
        Path::new(self.path).to_path_buf()
    }

    fn get_barcode_path(&self) -> std::path::PathBuf {
        Path::new(self.path).join("barcodes.csv")
    }

    fn get_historian(&self) -> Option<Rc<RefCell<Historian>>> {
        self.historian.clone()
    }

    fn write(
        &self,
        population: &Population<Store<S>>,
        generation: usize,
        compartment: usize,
    ) -> Result<String, std::io::Error> {
        let barcode = format!("sample_{generation}_{compartment}");

        // create output files
        let sample_path = Path::new(self.path).join(format!("{barcode}.fasta"));
        let mut samples_file = io::BufWriter::new(fs::File::create(sample_path).unwrap());

        // precompute sequences
        let sequences: HashMap<usize, Vec<u8>> = population
            .iter()
            .unique()
            .map(|haplotype_ref| {
                let sequence = haplotype_ref
                    .get_sequence()
                    .into_iter()
                    .map(|symbol| symbol.encode())
                    .collect();
                (haplotype_ref.get_id(), sequence)
            })
            .collect();

        // write to file
        for (haplotype_id, haplotype_ref) in population.iter().enumerate() {
            let head = format!(
                "sequence_id={};block_id={};compartment_id={};generation={}",
                haplotype_id,
                haplotype_ref.get_block_id(),
                compartment,
                generation
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
        Ok(barcode)
    }
}

impl<S: Symbol + 'static> SampleWriter<S> for CsvSampleWriter<'_> {
    fn get_experiment_name(&self) -> &str {
        self.simulation_name
    }

    fn get_path(&self) -> PathBuf {
        Path::new(self.path).to_path_buf()
    }

    fn get_barcode_path(&self) -> std::path::PathBuf {
        Path::new(self.path).join("barcodes.csv")
    }

    fn get_historian(&self) -> Option<Rc<RefCell<Historian>>> {
        self.historian.clone()
    }

    fn write(
        &self,
        population: &Population<Store<S>>,
        generation: usize,
        compartment: usize,
    ) -> Result<String, std::io::Error> {
        let barcode = format!("sample_{generation}_{compartment}");

        // create output file
        let sample_path = Path::new(self.path).join(format!("{barcode}.csv"));
        let mut samples_file = csv::WriterBuilder::new()
            .from_path(sample_path)
            .expect("Unable to open sample file.");

        // write header
        samples_file
            .write_record([
                "block_id",
                "compartment_id",
                "generation",
                "haplotype",
                "count",
            ])
            .expect("Unable to write header to samples file.");

        for (haplotype_ref, haplotype_count) in population.iter().counts() {
            samples_file
                .write_record(&[
                    haplotype_ref.get_block_id().to_string(),
                    compartment.to_string(),
                    generation.to_string(),
                    haplotype_ref.get_string(),
                    haplotype_count.to_string(),
                ])
                .expect("Unable to write to samples file.")
        }

        Ok(barcode)
    }
}

#[cfg(test)]
mod tests {
    use std::io::BufRead;

    use super::*;

    use crate::core::attributes::AttributeSetDefinition;
    use crate::core::haplotype::Wildtype;
    use crate::encoding::Nucleotide as Nt;
    use crate::references::HaplotypeRefTrait;

    fn get_population() -> Population<Store<Nt>> {
        let attribute_definition = AttributeSetDefinition::new();
        let wt1 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);

        population![&wt1; &wt1; &wt2; &wt2]
    }

    #[test]
    fn test_csv_write() {
        let tmp_dir = std::env::temp_dir();
        let path = tmp_dir.to_str().unwrap();
        let sample_writer = CsvSampleWriter::new::<Nt>("test_simulation", path, None).unwrap();

        let population = get_population();
        let barcode = sample_writer.write(&population, 0, 0).unwrap();

        let sample_path = Path::new(path).join(format!("{barcode}.csv"));

        assert!(sample_path.exists());

        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_path(sample_path)
            .unwrap();
        let first_record = reader.records().next().unwrap().unwrap();
        assert_eq!(first_record.len(), 5);
        assert_eq!(first_record[3], "wt".to_string());
        assert_eq!(first_record[4], "2".to_string());

        let second_record = reader.records().next().unwrap().unwrap();
        assert_eq!(second_record.len(), 5);
        assert_eq!(second_record[3], "wt".to_string());
        assert_eq!(second_record[4], "2".to_string());
    }

    #[test]
    fn test_fasta_write() {
        let tmp_dir = std::env::temp_dir();
        let path = tmp_dir.to_str().unwrap();
        let sample_writer = FastaSampleWriter::new::<Nt>("test_simulation", path, None).unwrap();

        let population = get_population();
        let barcode = sample_writer.write(&population, 0, 0).unwrap();

        let sample_path = Path::new(path).join(format!("{barcode}.fasta"));

        assert!(sample_path.exists());

        let file = fs::File::open(sample_path).unwrap();
        let reader = io::BufReader::new(file);
        let lines: Vec<String> = reader
            .lines()
            .map(|l| l.expect("Unable to read line"))
            .collect();

        assert_eq!(lines.len(), 8);
        assert_eq!(
            lines[0],
            format!(
                ">sequence_id=0;block_id={};compartment_id=0;generation=0",
                population.get(&0).get_block_id()
            )
        );
        assert_eq!(lines[1], "AAAAAAAAAA");
        assert_eq!(
            lines[2],
            format!(
                ">sequence_id=1;block_id={};compartment_id=0;generation=0",
                population.get(&1).get_block_id()
            )
        );
        assert_eq!(lines[3], "AAAAAAAAAA");
        assert_eq!(
            lines[4],
            format!(
                ">sequence_id=2;block_id={};compartment_id=0;generation=0",
                population.get(&2).get_block_id()
            )
        );
        assert_eq!(lines[5], "AAAAAAAAAA");
        assert_eq!(
            lines[6],
            format!(
                ">sequence_id=3;block_id={};compartment_id=0;generation=0",
                population.get(&3).get_block_id()
            )
        );
        assert_eq!(lines[7], "AAAAAAAAAA");
    }
}
