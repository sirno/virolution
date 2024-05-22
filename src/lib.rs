#![feature(let_chains)]
#![feature(get_mut_unchecked)]
#![feature(extract_if)]
#![feature(test)]

pub mod args;
pub mod barcode;
#[macro_use]
pub mod core;
pub mod config;
pub mod encoding;
pub mod readwrite;
pub mod references;
pub mod runner;
pub mod simulation;
pub mod stats;
