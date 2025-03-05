//! # Settings for the simulation binary
//!
//! Settings are used to configure the simulation when it is run as a program.
//! These settings are read from a YAML file and are used to initialize and run
//! the simulation. The settings are divided into two parts: parameters and
//! schedule. The parameters define the simulation parameters, such as the
//! mutation rate and the maximum virus population size. The schedule defines
//! events that occur throughout the simulation, such as transmission between
//! compartments, sampling of individuals, or environmental changes.
//!
//! ## Parameters
//!
//! The parameters define the simulation parameters and are given as a list of
//! [`Parameters`]. Some of the parameters require single values, while others
//! may require more complex structures. For example, the fitness model can
//! either be a single host model or a multi-host model. The corresponding MFED
//! may be a neutral model, a lognormal model, or a model based on a file.
//!
//! ```yaml
//! parameters:
//!  - mutation_rate: 1e-6
//!    recombination_rate: 0
//!    host_population_size: 10000
//!    basic_reproductive_number: 100.0
//!    max_population: 10000
//!    dilution: 0.02
//!    substitution_matrix:
//!      - [0.0, 1.0, 1.0, 1.0]
//!      - [1.0, 0.0, 1.0, 1.0]
//!      - [1.0, 1.0, 0.0, 1.0]
//!      - [1.0, 1.0, 1.0, 0.0]
//!    host_model: !SingleHost
//!      distribution: !Exponential
//!        weights:
//!          beneficial: 0.29
//!          deleterious: 0.51
//!          lethal: 0.2
//!          neutral: 0.0
//!        lambda_beneficial: 0.03
//!        lambda_deleterious: 0.21
//!      utility: !Algebraic
//!        upper: 1.5
//! ```
//!
//! All available options can be found by looking at the type definition of the
//! field. For example, one can find the available options for the
//! `host_model` field by following the type definition to the corresponding
//! enum: [`HostModel`]. The available options can then be defined in
//! YAML by using the corresponding tag, e.g., `!SingleHost` for a single host
//! model.
//!
//! ```yaml
//! host_model: !SingleHost
//!   distribution: !Neutral
//!   utility: !Linear
//! ```
//!
//! ## Schedule
//!
//! The schedule is defined as a list of events that occur throughout the
//! simulation. The events are defined by a generation expression, an event
//! name, and a value. The generation expression is either the generation where
//! the event should be executed, or in the case of recurring events an
//! expression that will be executed every time the expression evaluates to 0.
//! For example, the expression `{} % 200` will execute the event every 200
//! generations.
//!
//! ```yaml
//! schedule:
//!  - generation: "{} % 1"
//!    event: transmission
//!    value: "[[0.9, 0.1], [0.1, 0.9]]"
//!  - generation: "{} % 200"
//!    event: sample
//!    value: 1000
//! ```
//!

mod parameters;
mod schedule;
mod settings;

pub use parameters::{HostModel, HostFitness, Parameters};
pub use schedule::Schedule;
pub use settings::Settings;
