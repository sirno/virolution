use crate::haplotype::Haplotype;
use std::cell::RefCell;
use std::rc::{Rc, Weak};

#[derive(Clone)]
pub struct HaplotypeRef(pub Rc<RefCell<Haplotype>>);

impl std::ops::Deref for HaplotypeRef {
    // type Target = Rc<RefCell<Haplotype>>;
    type Target = Rc<RefCell<Haplotype>>;

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
        Self(Rc::new(RefCell::new(haplotype)))
        // Self(Arc::new(RwLock::new(haplotype)))
    }

    pub fn get_clone(&self) -> HaplotypeRef {
        HaplotypeRef(self.0.clone())
    }

    pub fn get_weak(&self) -> HaplotypeWeak {
        HaplotypeWeak(Rc::downgrade(&self.0))
        // HaplotypeWeak(Arc::downgrade(&self.0))
    }
}

#[derive(Clone)]
pub struct HaplotypeWeak(Weak<RefCell<Haplotype>>);

impl std::ops::Deref for HaplotypeWeak {
    type Target = Weak<RefCell<Haplotype>>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl HaplotypeWeak {
    pub fn upgrade(&self) -> Option<HaplotypeRef> {
        match self.0.upgrade() {
            Some(x) => Some(HaplotypeRef(x)),
            None => None,
        }
    }
}
