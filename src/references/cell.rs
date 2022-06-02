use crate::haplotype::Haplotype;
use derive_more::{Deref, DerefMut};
use std::cell::RefCell;
use std::rc::{Rc, Weak};

#[derive(Clone, Deref, DerefMut)]
pub struct HaplotypeRef(pub Rc<RefCell<Haplotype>>);

impl HaplotypeRef {
    pub fn new(haplotype: Haplotype) -> Self {
        Self(Rc::new(RefCell::new(haplotype)))
    }

    pub fn new_cyclic<F>(data_fn: F) -> Self
    where
        F: std::ops::Fn(&HaplotypeWeak) -> Haplotype,
    {
        Self(Rc::new_cyclic(|weak| {
            RefCell::new(data_fn(&HaplotypeWeak(weak.clone())))
        }))
    }

    #[inline]
    pub fn get_clone(&self) -> HaplotypeRef {
        HaplotypeRef(self.0.clone())
    }

    #[inline]
    pub fn get_weak(&self) -> HaplotypeWeak {
        HaplotypeWeak(Rc::downgrade(&self.0))
    }
}

impl PartialEq for HaplotypeRef {
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        Rc::ptr_eq(self, other)
    }
}

#[derive(Clone, Deref)]
pub struct HaplotypeWeak(Weak<RefCell<Haplotype>>);

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
    #[inline]
    fn eq(&self, other: &Self) -> bool {
        self.ptr_eq(other)
    }
}
