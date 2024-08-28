//! This module contains a fitness provider and any associated structures to load and initialize
//! fitness tables.

mod epistasis;
mod table;

pub mod init;
pub mod provider;
pub mod utility;

pub use provider::FitnessProvider;
pub use table::FitnessTable;
