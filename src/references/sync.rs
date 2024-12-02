#![allow(dead_code)]

use block_id::{Alphabet, BlockId};
use derive_more::{Deref, DerefMut};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Weak};

use crate::core::Haplotype;
use crate::encoding::Symbol;

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeRef<S: Symbol>(pub Arc<Haplotype<S>>);

thread_local! {
    pub static BLOCK_ID: BlockId<char> = BlockId::new(Alphabet::new(&("ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().collect::<Vec<char>>())), 0, 1);
}

impl<S: Symbol> HaplotypeRef<S> {
    pub fn new(haplotype: Haplotype<S>) -> Self {
        Self(Arc::new(haplotype))
    }

    pub fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::FnOnce(&HaplotypeWeak<S>) -> Haplotype<S>,
    {
        Self(Arc::new_cyclic(|weak| {
            data_fn(&HaplotypeWeak(weak.clone()))
        }))
    }

    #[inline]
    pub fn get_block_id(&self) -> String {
        let reference_ptr = Arc::as_ptr(&self.0) as u64;
        BLOCK_ID
            .with(|generator| generator.encode_string(reference_ptr))
            .expect("Unable to generate block ID")
    }

    #[inline]
    pub fn get_id(&self) -> usize {
        Arc::as_ptr(&self.0) as usize
    }

    #[inline]
    pub fn get_strong_count(&self) -> usize {
        Arc::strong_count(&self.0)
    }

    #[inline]
    pub fn get_weak(&self) -> HaplotypeWeak<S> {
        HaplotypeWeak(Arc::downgrade(&self.0))
    }

    #[inline]
    pub fn try_unwrap(&self) -> Option<Haplotype<S>> {
        Arc::try_unwrap(self.0.clone()).ok()
    }

    #[inline]
    pub fn as_ptr(&self) -> *const Haplotype<S> {
        Arc::as_ptr(&self.0)
    }
}

impl<S: Symbol> fmt::Debug for HaplotypeRef<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.get_string())
    }
}

impl<S: Symbol> PartialEq for HaplotypeRef<S> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(self, other)
    }
}

impl<S: Symbol> Eq for HaplotypeRef<S> {}

impl<S: Symbol> Hash for HaplotypeRef<S> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        Arc::as_ptr(&self.0).hash(state)
    }
}

impl<S: Symbol> Drop for HaplotypeRef<S> {
    fn drop(&mut self) {
        if Arc::strong_count(&self.0) == 2
            && let Some(mutant) = self.try_unwrap_mutant()
        {
            mutant.try_merge_node();
        }
    }
}

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeWeak<S: Symbol>(Weak<Haplotype<S>>);

impl<S: Symbol> HaplotypeWeak<S> {
    #[inline]
    pub fn upgrade(&self) -> Option<HaplotypeRef<S>> {
        self.0.upgrade().map(HaplotypeRef)
    }

    #[inline]
    pub fn exists(&self) -> bool {
        self.0.strong_count() > 0
    }

    #[inline]
    pub fn get_block_id(&self) -> String {
        let reference_ptr = Weak::as_ptr(&self.0) as u64;
        BLOCK_ID
            .with(|generator| generator.encode_string(reference_ptr))
            .expect("Unable to generate block ID")
    }

    #[inline]
    pub fn get_id(&self) -> usize {
        Weak::as_ptr(&self.0) as usize
    }
}

impl<S: Symbol> fmt::Debug for HaplotypeWeak<S> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.upgrade() {
            Some(reference) => write!(f, "{}", reference.get_string()),
            None => write!(f, "(None)"),
        }
    }
}

impl<S: Symbol> PartialEq for HaplotypeWeak<S> {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.ptr_eq(other)
    }
}
