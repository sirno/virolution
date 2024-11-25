use seq_io::fasta;
use seq_io::fasta::Record;

use crate::core::attributes::AttributeSet;
use crate::core::haplotype::{Haplotype, Wildtype};
use crate::encoding::Symbol;
use crate::errors::{Result, VirolutionError};
use crate::references::HaplotypeRef;

pub trait HaplotypeIO<S: Symbol> {
    /// Load wildtype from a fasta file
    fn load_wildtype_from_file(path: &str, attributes: AttributeSet<S>) -> Result<HaplotypeRef<S>>;

    /// Load wildtype from a sequence
    fn load_wildtype(sequence: Vec<S>, attributes: AttributeSet<S>) -> HaplotypeRef<S>;

    /// Read a sequence from a fasta file
    fn read_sequence_from_file(path: &str) -> Result<Vec<S>> {
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
            Some(Ok(sequence_record)) => Ok(decode_sequence::<S>(sequence_record)?),
        }
    }
}

fn decode_sequence<S: Symbol>(record: impl Record) -> Result<Vec<S>> {
    let sequence: Vec<S> = record
        .seq()
        .iter()
        .filter_map(|s| S::try_decode(s))
        .collect();
    if sequence.is_empty() {
        Err(VirolutionError::InitializationError(format!(
            "No valid symbols found in sequence: {}",
            record.id().unwrap_or("unknown")
        )))
    } else {
        Ok(sequence)
    }
}

impl<S: Symbol> HaplotypeIO<S> for Haplotype<S> {
    fn load_wildtype_from_file(path: &str, attributes: AttributeSet<S>) -> Result<HaplotypeRef<S>> {
        let sequence = Self::read_sequence_from_file(path)?;
        Ok(Wildtype::new(sequence, attributes))
    }

    fn load_wildtype(sequence: Vec<S>, attributes: AttributeSet<S>) -> HaplotypeRef<S> {
        Wildtype::new(sequence, attributes)
    }
}
