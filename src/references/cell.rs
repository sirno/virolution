#![allow(dead_code)]

use block_id::{Alphabet, BlockId};
use derive_more::{Deref, DerefMut};
use derive_where::derive_where;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::rc::{Rc, Weak};

use crate::core::Haplotype;
use crate::encoding::Symbol;

#[derive(Deref, DerefMut)]
#[derive_where(Clone)]
pub struct HaplotypeRef<S: Symbol>(pub Rc<Haplotype<S>>);

thread_local! {
    pub static BLOCK_ID: BlockId<char> = BlockId::new(Alphabet::new(&("ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().collect::<Vec<char>>())), 0, 1);
}

impl<S: Symbol> HaplotypeRef<S> {
    pub fn new(haplotype: Haplotype<S>) -> Self {
        Self(Rc::new(haplotype))
    }

    pub fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::FnOnce(&HaplotypeWeak<S>) -> Haplotype<S>,
    {
        Self(Rc::new_cyclic(|weak| data_fn(&HaplotypeWeak(weak.clone()))))
    }

    #[inline]
    pub fn get_strong_count(&self) -> usize {
        Rc::strong_count(&self.0)
    }

    #[inline]
    pub fn get_weak(&self) -> HaplotypeWeak<S> {
        HaplotypeWeak(Rc::downgrade(&self.0))
    }

    #[inline]
    pub fn get_block_id(&self) -> String {
        let reference_ptr = Rc::as_ptr(&self.0) as u64;
        BLOCK_ID.with(|generator| generator.encode_string(reference_ptr))
    }

    #[inline]
    pub fn get_id(&self) -> usize {
        Rc::as_ptr(&self.0) as usize
    }

    #[inline]
    pub fn try_unwrap(&self) -> Option<Haplotype<S>> {
        Rc::try_unwrap(self.0.clone()).ok()
    }

    #[inline]
    pub fn as_ptr(&self) -> *const Haplotype<S> {
        Rc::as_ptr(&self.0)
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
        Rc::ptr_eq(self, other)
    }
}

impl<S: Symbol> Eq for HaplotypeRef<S> {}

impl<S: Symbol> Hash for HaplotypeRef<S> {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        Rc::as_ptr(&self.0).hash(state);
    }
}

impl<S: Symbol> Drop for HaplotypeRef<S> {
    fn drop(&mut self) {
        if Rc::strong_count(&self.0) == 2
            && let Some(mutant) = self.try_unwrap_mutant()
        {
            mutant.try_merge_node();
        }
    }
}

#[derive(Deref, DerefMut)]
#[derive_where(Clone)]
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
        BLOCK_ID.with(|generator| generator.encode_string(reference_ptr))
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
