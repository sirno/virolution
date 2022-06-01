use crate::haplotype::Haplotype;
use derive_more::{Deref, DerefMut};
use std::sync::{Arc, Weak};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeRef(pub Arc<RwLock<Haplotype>>);

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
            if let Ok(guard) = self.0.try_write() {
                return guard;
            }
        }
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
}

impl PartialEq for HaplotypeWeak {
    fn eq(&self, other: &HaplotypeWeak) -> bool {
        self.ptr_eq(other)
    }
}
