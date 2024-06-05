//! This module contains reference implementations and chooses which implementation is visible
//! based on the feature flags.
mod cell;
mod desc;
mod sync;

#[cfg(feature = "parallel")]
pub use sync::{HaplotypeRef, HaplotypeWeak};

#[cfg(not(feature = "parallel"))]
pub use cell::{HaplotypeRef, HaplotypeWeak};

pub use desc::DescendantsCell;
