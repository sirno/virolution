use derive_more::Deref;
use parking_lot::Mutex;

use super::HaplotypeWeak;

#[derive(Debug, Deref)]
pub struct DescendantsCell(Mutex<Vec<HaplotypeWeak>>);

impl<'a> DescendantsCell {
    pub fn new() -> Self {
        Self(Mutex::new(Vec::new()))
    }

    pub fn from_iter(iter: impl IntoIterator<Item = HaplotypeWeak>) -> Self {
        Self(Mutex::new(iter.into_iter().collect()))
    }
}
