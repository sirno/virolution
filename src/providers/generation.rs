use std::sync::atomic::{AtomicUsize, Ordering};

use crate::core::{AttributeProvider, AttributeProviderType, AttributeValue};
use crate::encoding::Symbol;
use crate::errors::Result;
use crate::references::HaplotypeRef;

#[derive(Debug)]
pub struct Generation {
    generation: AtomicUsize,
}

/// A provider for a generation attribute.
///
/// The `Generation` provider is intended to track discrete time evolution in a simulation.
impl Generation {
    pub fn new(generation: usize) -> Self {
        Self {
            generation: AtomicUsize::new(generation),
        }
    }

    pub fn get(&self) -> usize {
        self.generation.load(Ordering::Relaxed)
    }

    pub fn set(&self, generation: usize) {
        self.generation.store(generation, Ordering::Relaxed);
    }

    pub fn increment(&self) {
        self.generation.fetch_add(1, Ordering::Relaxed);
    }
}

impl<S: Symbol> AttributeProvider<S> for Generation {
    fn name(&self) -> &'static str {
        "generation"
    }

    fn compute(&self, _haplotype: &Option<HaplotypeRef<S>>) -> AttributeValue {
        AttributeValue::U64(self.generation.load(Ordering::Relaxed) as u64)
    }

    fn get_provider_type(&self) -> AttributeProviderType {
        AttributeProviderType::Eager
    }

    fn write(&self, _path: &std::path::Path) -> Result<()> {
        unimplemented!()
    }
}
