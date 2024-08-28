//! This module contains the core datatypes of the library.

mod ancestry;
mod historian;

#[macro_use]
pub mod population;
pub mod fitness;
pub mod haplotype;

pub use ancestry::Ancestry;
pub use fitness::FitnessProvider;
pub use haplotype::Haplotype;
pub use historian::Historian;
pub use population::Population;
