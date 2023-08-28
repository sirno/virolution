#![feature(let_chains)]
#![feature(get_mut_unchecked)]
#![feature(extract_if)]
#![feature(test)]

pub mod args;
pub mod barcode;
pub mod fitness;
pub mod haplotype;
pub mod historian;
#[macro_use]
pub mod population;
pub mod ancestry;
pub mod references;
pub mod sample_writer;
pub mod simulation;
pub mod simulation_parameters;
pub mod simulation_plan;
pub mod simulation_settings;
