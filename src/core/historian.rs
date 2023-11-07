//! Historian --- stores events in memory to keep elements from being dropped
//!
//! The historian is a data structure that stores events in memory. Its primary
//! purpose is to keep elements from being dropped if only weak references to
//! them exist. This is useful for keeping track of the history of a simulation
//! without having to keep the entire history in memory of unobserved elements.
//!

use crate::core::Population;
use std::fmt;

#[derive(Debug)]
pub struct Historian {
    history: Vec<Box<dyn HistoricalEvent>>,
}

trait HistoricalEvent: std::fmt::Debug {}

#[derive(Debug)]
struct SampleEvent {
    generation: usize,
    compartment: usize,
    sample: Population,
}

impl HistoricalEvent for SampleEvent {}

impl fmt::Display for SampleEvent {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "SampleEvent(generation={},compartment={},length={})",
            self.generation,
            self.compartment,
            self.sample.len()
        )
    }
}

impl SampleEvent {
    fn new(generation: usize, compartment: usize, sample: Population) -> Self {
        Self {
            generation,
            compartment,
            sample,
        }
    }
}

impl Historian {
    pub fn new() -> Self {
        Self { history: vec![] }
    }

    pub fn record_sample(&mut self, generation: usize, compartment: usize, sample: Population) {
        self.history
            .push(Box::new(SampleEvent::new(generation, compartment, sample)));
    }
}

impl Default for Historian {
    fn default() -> Self {
        Self::new()
    }
}
