use crate::core::haplotype::{Haplotype, Symbol, Wildtype};
use crate::encoding::STRICT_DECODE;
use crate::references::HaplotypeRef;
use seq_io::fasta;
use seq_io::fasta::Record;
use std::panic::catch_unwind;

pub trait HaplotypeIO {
    fn load_wildtype(path: &str) -> Result<HaplotypeRef, fasta::Error>;
}

impl HaplotypeIO for Haplotype {
    fn load_wildtype(path: &str) -> Result<HaplotypeRef, fasta::Error> {
        let mut reader = fasta::Reader::from_path(path)?;
        let sequence: Vec<Symbol> = reader
            .next()
            .unwrap()
            .expect("Unable to read sequence.")
            .seq()
            .iter()
            .filter(|&&enc| enc != 0x0au8)
            .map(|enc| match catch_unwind(|| Some(STRICT_DECODE[enc])) {
                Ok(result) => result,
                Err(_) => panic!("Unable to decode literal {enc}."),
            })
            .collect();
        Ok(Wildtype::new(sequence))
    }
}
