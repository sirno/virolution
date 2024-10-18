use crate::encoding::Symbol;
use derive_more::Deref;
use parking_lot::Mutex;

use super::HaplotypeWeak;

#[derive(Debug, Deref)]
pub struct DescendantsCell<S: Symbol>(Mutex<Vec<HaplotypeWeak<S>>>);

impl<S: Symbol> DescendantsCell<S> {
    pub fn new() -> Self {
        Self(Mutex::new(Vec::new()))
    }
}

impl<S: Symbol> FromIterator<HaplotypeWeak<S>> for DescendantsCell<S> {
    fn from_iter<I: IntoIterator<Item = HaplotypeWeak<S>>>(iter: I) -> Self {
        Self(Mutex::new(iter.into_iter().collect()))
    }
}

impl<S: Symbol> Default for DescendantsCell<S> {
    fn default() -> Self {
        Self::new()
    }
}
