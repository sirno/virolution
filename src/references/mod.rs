mod cell;
mod sync;

#[cfg(feature = "parallel")]
pub use sync::{HaplotypeRef, HaplotypeWeak};

#[cfg(not(feature = "parallel"))]
pub use cell::{HaplotypeRef, HaplotypeWeak};
