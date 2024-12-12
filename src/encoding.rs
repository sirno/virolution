//! Encoding and decoding definitions for Symbol implementations.

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Nucleotide {
    A,
    T,
    C,
    G,
}

pub trait Symbol:
    std::marker::Sized
    + Copy
    + Clone
    + Send
    + Sync
    + std::fmt::Debug
    + std::cmp::PartialEq
    + std::cmp::Eq
    + std::hash::Hash
    + std::fmt::Display
    + 'static
{
    const SIZE: usize;
    fn decode(s: &u8) -> Self;
    fn try_decode(s: &u8) -> Option<Self>;
    fn encode(&self) -> u8;
    fn index(&self) -> &usize;
}

impl Symbol for Nucleotide {
    const SIZE: usize = 4;

    fn decode(s: &u8) -> Self {
        match s {
            // ATCG | atcg | 0123 -> Nucleotide
            0x41 | 0x61 | 0x30 | 0x00 => Nucleotide::A,
            0x54 | 0x74 | 0x31 | 0x01 => Nucleotide::T,
            0x43 | 0x63 | 0x32 | 0x02 => Nucleotide::C,
            0x47 | 0x67 | 0x33 | 0x03 => Nucleotide::G,
            _ => panic!("Invalid nucleotide symbol: {}", s),
        }
    }

    fn try_decode(s: &u8) -> Option<Self> {
        match s {
            // ATCG | atcg | 0123 -> Nucleotide
            0x41 | 0x61 | 0x30 | 0x00 => Some(Nucleotide::A),
            0x54 | 0x74 | 0x31 | 0x01 => Some(Nucleotide::T),
            0x43 | 0x63 | 0x32 | 0x02 => Some(Nucleotide::C),
            0x47 | 0x67 | 0x33 | 0x03 => Some(Nucleotide::G),
            _ => None,
        }
    }

    fn encode(&self) -> u8 {
        match self {
            Nucleotide::A => 0x41,
            Nucleotide::T => 0x54,
            Nucleotide::C => 0x43,
            Nucleotide::G => 0x47,
        }
    }

    fn index(&self) -> &usize {
        match self {
            Nucleotide::A => &0,
            Nucleotide::T => &1,
            Nucleotide::C => &2,
            Nucleotide::G => &3,
        }
    }
}

impl std::fmt::Display for Nucleotide {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "{}", self.encode() as char)
    }
}
