use derive_more::Deref;
use parking_lot::Mutex;

use super::HaplotypeWeak;

#[derive(Debug, Deref)]
pub struct DescendantsCell(Mutex<Vec<HaplotypeWeak>>);

impl DescendantsCell {
    pub fn new() -> Self {
        Self(Mutex::new(Vec::new()))
    }
}

impl FromIterator<HaplotypeWeak> for DescendantsCell {
    fn from_iter<I: IntoIterator<Item=HaplotypeWeak>>(iter: I) -> Self {
        Self(Mutex::new(iter.into_iter().collect()))
    }
}

impl Default for DescendantsCell {
    fn default() -> Self {
        Self::new()
    }
}
