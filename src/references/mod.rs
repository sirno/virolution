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

pub trait HaplotypeRefTrait:
    Clone
    + std::ops::Deref<Target = crate::core::haplotype::Haplotype<Self::Symbol>>
    + std::fmt::Debug
    + std::cmp::PartialEq
    + std::cmp::Eq
    + std::hash::Hash
{
    type Symbol: crate::encoding::Symbol;
    type Weak: HaplotypeWeakTrait<Symbol = Self::Symbol>;

    fn new(haplotype: crate::core::haplotype::Haplotype<Self::Symbol>) -> Self;
    fn new_cyclic<
        F: std::ops::FnOnce(&Self::Weak) -> crate::core::haplotype::Haplotype<Self::Symbol>,
    >(
        data_fn: F,
    ) -> Self;
    fn get_strong_count(&self) -> usize;
    fn get_weak(&self) -> Self::Weak;
    fn get_block_id(&self) -> String;
    fn get_id(&self) -> usize;
    fn try_unwrap(&self) -> Option<crate::core::haplotype::Haplotype<Self::Symbol>>;
    fn as_ptr(&self) -> *const crate::core::haplotype::Haplotype<Self::Symbol>;
}

pub trait HaplotypeWeakTrait: Clone + Sized {
    type Symbol: crate::encoding::Symbol;
    type Ref: HaplotypeRefTrait<Symbol = Self::Symbol>;

    fn upgrade(&self) -> Option<Self::Ref>;
    fn exists(&self) -> bool;
    fn get_block_id(&self) -> String;
    fn get_id(&self) -> usize;
    fn as_ptr(&self) -> *const crate::core::haplotype::Haplotype<Self::Symbol>;
}
