//! This module contains a fitness provider and any associated structures to load and initialize
//! fitness tables.

pub(crate) mod epistasis;
pub(crate) mod table;
pub(crate) mod utility;

pub use epistasis::EpistasisTable;
pub use table::FitnessTable;
pub use utility::UtilityFunction;
