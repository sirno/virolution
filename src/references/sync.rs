use crate::haplotype::Haplotype;
use block_id::{Alphabet, BlockId};
use derive_more::{Deref, DerefMut};
use std::fmt;
use std::hint;
use std::sync::{Arc, Weak};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeRef(pub Arc<RwLock<Haplotype>>);

thread_local! {
    pub static BLOCK_ID: BlockId<char> = BlockId::new(Alphabet::new(&("ABCDEFGHIJKLMNOPQRSTUVWXYZ".chars().collect::<Vec<char>>())), 0, 1);
}

fn print_reference_weak(
    reference_weak: &HaplotypeWeak,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    match reference_weak.upgrade() {
        Some(reference) => write!(formatter, "{}", reference.borrow().get_string()),
        None => write!(formatter, "None"),
    }
}

impl HaplotypeRef {
    pub fn new(haplotype: Haplotype) -> Self {
        Self(Arc::new(RwLock::new(haplotype)))
    }

    pub fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::Fn(&HaplotypeWeak) -> Haplotype,
    {
        Self(Arc::new_cyclic(|weak| {
            RwLock::new(data_fn(&HaplotypeWeak(weak.clone())))
        }))
    }

    #[inline]
    pub fn get_id(&self) -> String {
        let reference_ptr = Arc::as_ptr(&self.0) as u64;
        BLOCK_ID.with(|generator| generator.encode_string(reference_ptr))
    }

    #[inline]
    pub fn get_clone(&self) -> HaplotypeRef {
        HaplotypeRef(self.0.clone())
    }

    #[inline]
    pub fn get_weak(&self) -> HaplotypeWeak {
        HaplotypeWeak(Arc::downgrade(&self.0))
    }

    #[inline]
    pub fn borrow(&self) -> RwLockReadGuard<'_, Haplotype> {
        self.0.read().unwrap()
    }

    #[inline]
    pub fn borrow_mut(&self) -> RwLockWriteGuard<'_, Haplotype> {
        loop {
            hint::spin_loop();
            if let Ok(guard) = self.0.try_write() {
                return guard;
            }
        }
    }
}

impl fmt::Debug for HaplotypeRef {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.borrow().get_string())
    }
}

impl PartialEq for HaplotypeRef {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        Arc::ptr_eq(self, other)
    }
}

#[derive(Clone, Deref)]
pub struct HaplotypeWeak(Weak<RwLock<Haplotype>>);

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
    pub fn get_id(&self) -> String {
        let reference_ptr = Weak::as_ptr(&self.0) as u64;
        BLOCK_ID.with(|generator| generator.encode_string(reference_ptr))
    }
}

impl fmt::Debug for HaplotypeWeak {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self.upgrade() {
            Some(reference) => write!(f, "{}", reference.borrow().get_string()),
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
