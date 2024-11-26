//! IO traits for reading and writing haplotypes and populations.

mod haplotype;
mod population;
mod sample_writer;

pub use haplotype::HaplotypeIO;
pub use population::PopulationIO;
pub use sample_writer::{CsvSampleWriter, FastaSampleWriter, SampleWriter};
