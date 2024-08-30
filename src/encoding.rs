//!
//! This module contains the encoding and decoding definitions for nucleotide
//! symbols.
//!

use phf::phf_map;

pub static STRICT_ENCODE: phf::Map<u8, u8> = phf_map! {
    0x00u8 => 0x41,
    0x01u8 => 0x54,
    0x02u8 => 0x43,
    0x03u8 => 0x47,
};

pub static STRICT_DECODE: phf::Map<u8, u8> = phf_map! {
    0x41u8 => 0x00,
    0x54u8 => 0x01,
    0x43u8 => 0x02,
    0x47u8 => 0x03,
};

pub static DECODE: phf::Map<u8, u8> = phf_map! {
    // ATCG -> 0123
    0x41u8 => 0x00,
    0x54u8 => 0x01,
    0x43u8 => 0x02,
    0x47u8 => 0x03,
    // atcg -> 0123
    0x61u8 => 0x00,
    0x74u8 => 0x01,
    0x63u8 => 0x02,
    0x67u8 => 0x03,
    // 0123 -> 0123
    0x30u8 => 0x00,
    0x31u8 => 0x01,
    0x32u8 => 0x02,
    0x33u8 => 0x03,
};