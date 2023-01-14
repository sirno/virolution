#![allow(dead_code)]

use crate::haplotype::Haplotype;
use block_id::{Alphabet, BlockId};
use derive_more::{Deref, DerefMut};
use std::fmt;
use std::hash::{Hash, Hasher};
use std::rc::{Rc, Weak};

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeRef(pub Rc<Haplotype>);

thread_local! {
    pub static BLOCK_ID: BlockId<char> = BlockId::new(Alphabet::new(&("ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().collect::<Vec<char>>())), 0, 1);
}

impl HaplotypeRef {
    pub fn new(haplotype: Haplotype) -> Self {
        Self(Rc::new(haplotype))
    }

    pub fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::Fn(&HaplotypeWeak) -> Haplotype,
    {
        Self(Rc::new_cyclic(|weak| data_fn(&HaplotypeWeak(weak.clone()))))
    }

    #[inline]
    pub fn get_clone(&self) -> HaplotypeRef {
        HaplotypeRef(self.0.clone())
    }

    #[inline]
    pub fn get_weak(&self) -> HaplotypeWeak {
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
}

impl fmt::Debug for HaplotypeRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.get_string())
    }
}

impl PartialEq for HaplotypeRef {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        Rc::ptr_eq(self, other)
    }
}

impl Eq for HaplotypeRef {}

impl Hash for HaplotypeRef {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        Rc::as_ptr(&self.0).hash(state);
    }
}

#[derive(Clone, Deref)]
pub struct HaplotypeWeak(Weak<Haplotype>);

impl HaplotypeWeak {
    #[inline]
    pub fn upgrade(&self) -> Option<HaplotypeRef> {
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

impl fmt::Debug for HaplotypeWeak {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.upgrade() {
            Some(reference) => write!(f, "{}", reference.get_string()),
            None => write!(f, "(None)"),
        }
    }
}

impl PartialEq for HaplotypeWeak {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.ptr_eq(other)
    }
}
