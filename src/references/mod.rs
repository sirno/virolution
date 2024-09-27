//! This module contains reference implementations and chooses which implementation is visible
//! based on the feature flags.
mod cell;
mod desc;
mod sync;

pub use sync::{HaplotypeRef, HaplotypeWeak};

pub use desc::DescendantsCell;
