//! Configuration data structures for simulation setups.

mod parameters;
mod schedule;
mod settings;

pub use parameters::{FitnessModelField, Parameters};
pub use schedule::Schedule;
pub use settings::{Settings, SettingsError};
