# virolution

[![Crates.io Version](https://img.shields.io/crates/v/virolution?style=for-the-badge)](https://crates.io/crates/virolution)
[![docs.rs](https://img.shields.io/docsrs/virolution?style=for-the-badge)](https://docs.rs/virolution)
[![GitHub Branch](https://img.shields.io/github/check-runs/sirno/virolution/main?style=for-the-badge)](https://github.com/sirno/virolution/actions)
<!-- cargo-rdme start -->

## Virolution - Agent-Based Simulation of Viral Evolution

Virolution is a simulation framework that is powered by a sparse genetic representation of
viral genomes. It is designed to simulate the evolution of large viral populations with large
genomes. The simulations core data structures are thread-safe and generic to allow for
adaptation to a wide range of viral or other evolutionary systems.

### How to use Virolution

There are two ways of using Virolution:

- As a binary: You can use the binary to run simulations with a given configuration file. This
  will run a simulation that approximates the evolution of a viral population in discrete
  generations in potentially multiple compartments. The simulation can be configured to use
  different fitness landscapes, mutation rates, and population sizes.

- As a library: You can use the core data structures to build your own simulations. This allows
  you to create custom simulations that are not possible with the binary. Using the libraries
  traits you can implement your own fitness landscapes, attribute providers, and host systems.
  You can also use non-discrete time steps and implement your own selection mechanisms.

### Installation

To install the binary, you can use the following command:

```shell
cargo install virolution
cargo install virolution --features parallel
```

To use the library, you can add the following to your `Cargo.toml`:

```toml
[dependencies]
virolution = "0.5"
```

### Basic Usage as a Binary

To find all available options and their default values, you can run:

```shell
virolution --help
```

The simulation can be run with a sequence file and a configuration file. A configuration file
provides the parameters for the simulation environment and the schedule for any events during
the simulation. It may look like this:

```yaml
parameters:
 - mutation_rate: 1e-6
   recombination_rate: 0
   host_population_size: 10000
   basic_reproductive_number: 100.0
   max_population: 10000
   dilution: 0.02
   substitution_matrix:
     - [0.0, 1.0, 1.0, 1.0]
     - [1.0, 0.0, 1.0, 1.0]
     - [1.0, 1.0, 0.0, 1.0]
     - [1.0, 1.0, 1.0, 0.0]
   host_model: !SingleHost
     distribution: !Exponential
       weights:
         beneficial: 0.29
         deleterious: 0.51
         lethal: 0.2
         neutral: 0.0
       lambda_beneficial: 0.03
       lambda_deleterious: 0.21
     utility: !Algebraic
       upper: 1.5
schedule:
 - generation: "{} % 1"
   event: transmission
   value: "[[0.9, 0.1], [0.1, 0.9]]"
 - generation: "{} % 200"
   event: sample
   value: 1000
```

This configuration file will run a simulation with a mutation rate of 1e-6, a host population
size of 10000, and a basic reproductive number of 100.0. The substitution matrix is a 4x4
matrix with 0.0 on the diagonal and 1.0 elsewhere, meaning that all mutations are equally
likely. The fitness model is a single host model with an exponential distribution of fitness
effects, with a utility function that caps the fitness of the otherwise multiplicative fitness
at 1.5.

In the schedule, we define two events: transmission and sample. The transmission event will
occur every generation and will use a 2x2 matrix to determine the probability of transmission
between two compartments. The sample event will occur every 200 generations and will sample
1000 individuals from the population. The `{}` in the generation field will be replaced with
the current generation number and an event will be executed whenever the expression evaluates
to 0.

There are more parameters and events that can be configured in the schedule. For more
information on the configuration, see the
[documentation](https://docs.rs/virolution/latest/virolution/config/index.html).

When the simulation is complete (reached its maximum generation), the final population will be
stored in a file called `final.<i>.csv`, where `<i>` is the compartment number. The file will
contain all haplotypes and their respective counts.

### Basic Usage as a Library

To use Virolution as a library, you can simply create your own binary in Rust and use the
implementations provided by the library.

A wildtype haplotype is created from a sequence and an attribute definition. From it, mutants
and recombinants can be derived:

```rust
use std::sync::Arc;

use virolution::core::*;
use virolution::init::fitness::*;
use virolution::core::fitness::UtilityFunction;
use virolution::providers::FitnessProvider;
use virolution::encoding::Nucleotide as Nt;

let sequence = vec![Nt::A; 4];
let mut attribute_definitions = AttributeSetDefinition::new();
attribute_definitions.register(Arc::new(
    FitnessProvider::from_model(
        "fitness",
        &sequence,
        &FitnessModel::new(FitnessDistribution::Neutral, UtilityFunction::Linear),
    )
    .expect("Failed to create fitness table"),
));

let wildtype = Haplotype::new(sequence, &attribute_definitions);
let mutant1 = wildtype.create_descendant(vec![2], vec![Nt::T]);
let mutant2 = wildtype.create_descendant(vec![1, 2], vec![Nt::G, Nt::C]);
let recombinant = Haplotype::create_recombinant(&mutant1, &mutant2, 0, 2);

assert_eq!(wildtype.get_sequence(), vec![Nt::A, Nt::A, Nt::A, Nt::A]);
assert_eq!(mutant1.get_sequence(), vec![Nt::A, Nt::A, Nt::T, Nt::A]);
assert_eq!(mutant2.get_sequence(), vec![Nt::A, Nt::G, Nt::C, Nt::A]);
assert_eq!(recombinant.get_sequence(), vec![Nt::A, Nt::A, Nt::C, Nt::A]);
```

<!-- cargo-rdme end -->
