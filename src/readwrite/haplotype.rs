use crate::core::haplotype::{Haplotype, Symbol, Wildtype};
use crate::encoding::STRICT_DECODE;
use crate::references::HaplotypeRef;
use seq_io::fasta;
use seq_io::fasta::Record;

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
            .filter(|&&enc| enc != 0x0au8 && enc != 0x0du8)
            .map(|enc| match STRICT_DECODE.get(enc) {
                Some(result) => Some(*result),
                None => panic!("Unable to decode literal {enc}."),
            })
            .collect();
        Ok(Wildtype::new(sequence))
    }
}
