use seq_io::fasta;
use seq_io::fasta::Record;

use crate::core::haplotype::{Haplotype, Wildtype};
use crate::encoding::Symbol;
use crate::errors::{Result, VirolutionError};
use crate::references::HaplotypeRef;

pub trait HaplotypeIO<S: Symbol> {
    fn load_wildtype(path: &str) -> Result<HaplotypeRef<S>>;
}

impl<S: Symbol> HaplotypeIO<S> for Haplotype<S> {
    fn load_wildtype(path: &str) -> Result<HaplotypeRef<S>> {
        let mut reader = fasta::Reader::from_path(path).map_err(|_| {
            VirolutionError::InitializationError(format!(
                "Unable create file reader for fasta file: {path}"
            ))
        })?;

        match reader.next() {
            None => Err(VirolutionError::InitializationError(format!(
                "No sequence found in fasta file: {path}"
            ))),
            Some(Err(_)) => Err(VirolutionError::InitializationError(format!(
                "Unable to read sequence from fasta file: {path}"
            ))),
            Some(Ok(sequence_record)) => {
                let sequence: Vec<S> = sequence_record
                    .seq()
                    .iter()
                    .filter_map(|s| S::try_decode(s))
                    .collect();
                Ok(Wildtype::new(sequence))
            }
        }
    }
}
