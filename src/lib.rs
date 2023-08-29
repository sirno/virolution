#![feature(let_chains)]
#![feature(get_mut_unchecked)]
#![feature(extract_if)]
#![feature(test)]

pub mod args;
pub mod barcode;
#[macro_use]
pub mod core;
pub mod config;
pub mod references;
pub mod runner;
pub mod sample_writer;
pub mod simulation;
