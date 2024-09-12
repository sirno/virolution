use seq_io::fasta;
use seq_io::fasta::Record;

use crate::core::haplotype::{Haplotype, Symbol, Wildtype};
use crate::encoding::STRICT_DECODE;
use crate::errors::{Result, VirolutionError};
use crate::references::HaplotypeRef;

pub trait HaplotypeIO {
    fn load_wildtype(path: &str) -> Result<HaplotypeRef>;
}

impl HaplotypeIO for Haplotype {
    fn load_wildtype(path: &str) -> Result<HaplotypeRef> {
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
                let sequence: Vec<Symbol> = sequence_record
                    .seq()
                    .iter()
                    .filter_map(|enc| STRICT_DECODE.get(enc).copied().map(Some))
                    .collect();
                Ok(Wildtype::new(sequence))
            }
        }
    }
}
