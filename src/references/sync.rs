use crate::haplotype::Haplotype;
use std::sync::{Arc, Weak};
use std::sync::{RwLock, RwLockReadGuard, RwLockWriteGuard};

#[derive(Clone)]
pub struct HaplotypeRef(pub Arc<RwLock<Haplotype>>);

impl std::ops::Deref for HaplotypeRef {
    type Target = Arc<RwLock<Haplotype>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl std::ops::DerefMut for HaplotypeRef {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
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

    pub fn get_clone(&self) -> HaplotypeRef {
        HaplotypeRef(self.0.clone())
    }

    pub fn get_weak(&self) -> HaplotypeWeak {
        HaplotypeWeak(Arc::downgrade(&self.0))
    }

    pub fn borrow(&self) -> RwLockReadGuard<'_, Haplotype> {
        self.0.read().unwrap()
    }

    pub fn borrow_mut(&self) -> RwLockWriteGuard<'_, Haplotype> {
        loop {
            if let Ok(guard) = self.0.try_write() {
                return guard;
            }
        }
    }
}

#[derive(Clone)]
pub struct HaplotypeWeak(Weak<RwLock<Haplotype>>);

impl std::ops::Deref for HaplotypeWeak {
    type Target = Weak<RwLock<Haplotype>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl HaplotypeWeak {
    pub fn upgrade(&self) -> Option<HaplotypeRef> {
        self.0.upgrade().map(HaplotypeRef)
    }
}
