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

pub trait Identifiable {
    fn get_id(&self) -> usize;
}

pub trait DerefHaplotype<S: crate::encoding::Symbol> {
    fn deref_haplotype(&self) -> &crate::core::Haplotype<S>;
}

pub trait HaplotypeTraitBound<S: crate::encoding::Symbol>:
    Clone
    + Identifiable
    + DerefHaplotype<S>
    + std::fmt::Debug
    + std::cmp::PartialEq
    + std::cmp::Eq
    + std::hash::Hash
{
}
