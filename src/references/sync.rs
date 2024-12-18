#![allow(dead_code)]

use block_id::{Alphabet, BlockId};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, Weak};

use crate::core::Haplotype;
use crate::encoding::Symbol;

use super::{HaplotypeRefTrait, HaplotypeWeakTrait};

pub struct HaplotypeRef<S: Symbol>(pub Arc<Haplotype<S>>);

thread_local! {
    pub static BLOCK_ID: BlockId<char> = BlockId::new(Alphabet::new(&("ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().collect::<Vec<char>>())), 0, 1);
}

impl<S: Symbol> HaplotypeRefTrait for HaplotypeRef<S> {
    type Symbol = S;
    type Weak = HaplotypeWeak<Self::Symbol>;

    fn new(haplotype: Haplotype<S>) -> Self {
        Self(Arc::new(haplotype))
    }

    fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::FnOnce(&HaplotypeWeak<S>) -> Haplotype<S>,
    {
        Self(Arc::new_cyclic(|weak| {
            data_fn(&HaplotypeWeak(weak.clone()))
        }))
    }

    #[inline]
    fn get_block_id(&self) -> String {
        let reference_ptr = Arc::as_ptr(&self.0) as u64;
        BLOCK_ID
            .with(|generator| generator.encode_string(reference_ptr))
            .expect("Unable to generate block ID")
    }

    #[inline]
    fn get_id(&self) -> usize {
        Arc::as_ptr(&self.0) as usize
    }

    #[inline]
    fn get_strong_count(&self) -> usize {
        Arc::strong_count(&self.0)
    }

    #[inline]
    fn get_weak(&self) -> HaplotypeWeak<S> {
        HaplotypeWeak(Arc::downgrade(&self.0))
    }

    #[inline]
    fn try_unwrap(&self) -> Option<Haplotype<S>> {
        Arc::try_unwrap(self.0.clone()).ok()
    }

    #[inline]
    fn as_ptr(&self) -> *const Haplotype<S> {
        Arc::as_ptr(&self.0)
    }
}

impl<S: Symbol> std::ops::Deref for HaplotypeRef<S> {
    type Target = Haplotype<S>;

    fn deref(&self) -> &Self::Target {
        &self.0
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
        Arc::ptr_eq(&self.0, &other.0)
    }
}

impl<S: Symbol> Eq for HaplotypeRef<S> {}

impl<S: Symbol> Hash for HaplotypeRef<S> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        Arc::as_ptr(&self.0).hash(state)
    }
}

impl<S: Symbol> Clone for HaplotypeRef<S> {
    #[inline]
    fn clone(&self) -> Self {
        Self(self.0.clone())
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

pub struct HaplotypeWeak<S: Symbol>(Weak<Haplotype<S>>);

impl<S: Symbol> HaplotypeWeakTrait for HaplotypeWeak<S> {
    type Symbol = S;
    type Ref = HaplotypeRef<Self::Symbol>;

    #[inline]
    fn upgrade(&self) -> Option<HaplotypeRef<S>> {
        self.0.upgrade().map(HaplotypeRef)
    }

    #[inline]
    fn exists(&self) -> bool {
        self.0.strong_count() > 0
    }

    #[inline]
    fn get_block_id(&self) -> String {
        let reference_ptr = Weak::as_ptr(&self.0) as u64;
        BLOCK_ID
            .with(|generator| generator.encode_string(reference_ptr))
            .expect("Unable to generate block ID")
    }

    #[inline]
    fn get_id(&self) -> usize {
        Weak::as_ptr(&self.0) as usize
    }

    #[inline]
    fn as_ptr(&self) -> *const Haplotype<S> {
        Weak::as_ptr(&self.0)
    }
}

impl<S: Symbol> Clone for HaplotypeWeak<S> {
    fn clone(&self) -> Self {
        Self(self.0.clone())
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
        self.0.ptr_eq(&other.0)
    }
}
