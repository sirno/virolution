use crate::haplotype::Haplotype;
use derive_more::{Deref, DerefMut};
use std::cell::RefCell;
use std::rc::{Rc, Weak};

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeRef(pub Rc<RefCell<Haplotype>>);

impl HaplotypeRef {
    pub fn new(haplotype: Haplotype) -> Self {
        Self(Rc::new(RefCell::new(haplotype)))
        // Self(Arc::new(RwLock::new(haplotype)))
    }

    pub fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::Fn(&HaplotypeWeak) -> Haplotype,
    {
        Self(Rc::new_cyclic(|weak| {
            RefCell::new(data_fn(&HaplotypeWeak(weak.clone())))
        }))
    }

    pub fn get_clone(&self) -> HaplotypeRef {
        HaplotypeRef(self.0.clone())
    }

    pub fn get_weak(&self) -> HaplotypeWeak {
        HaplotypeWeak(Rc::downgrade(&self.0))
    }
}

#[derive(Clone, Deref)]
pub struct HaplotypeWeak(Weak<RefCell<Haplotype>>);

impl HaplotypeWeak {
    pub fn upgrade(&self) -> Option<HaplotypeRef> {
        self.0.upgrade().map(HaplotypeRef)
    }
}
