# Virolution

This tool simulates a virus population under static selective pressure with
single-nucleotide polymorphism and recombination in discrete generations.  The
host population is kept constant for each generation, while the virus population
may change in number and composition. Multiple simulations can be run in
parallel with controlled mixing. As input, the simulation requires model
parameters, a wildtype sequence and a file with a schedule for all migration and
sampling events. As output, multiple fasta files with samples of the populations
at each sampling event are generated.

## Details

In each generation every virus can infect a randomly selected host with
probability $r$. If multiple viruses infect the same host then each pair of
infectants will recombine with probability $p_r$. Next each of the infectants
can mutate on each site with probability $p_m$. Lastly, each of the infectants
will generate offspring. The number of offspring it can generate is drawn from
a Poisson-distribution, where the expectation is the product of fitness and
basic reproductive number, divided by the number of infectants in their host.

## Simulation Parameters

```rust
pub struct SimulationParameters {
    pub mutation_rate: f64,
    pub recombination_rate: f64,
    pub host_population_size: usize,
    pub infection_fraction: f64,
    pub basic_reproductive_number: f64,
    pub max_population: usize,
    pub dilution: f64,
    pub substitution_matrix: [[f64; 4]; 4],
    pub fitness_model: FitnessModel,
}
```

### Example

```yaml
---
simulation_parameters:
  mutation_rate: 1e-6
  recombination_rate: 1e-4
  host_population_size: 100000000
  infection_fraction: 0.7
  basic_reproductive_number: 100.0
  max_population: 100000000
  dilution: 0.02
  substitution_matrix:
    - [0.0, 1.0, 1.0, 1.0]
    - [1.0, 0.0, 1.0, 1.0]
    - [1.0, 1.0, 0.0, 1.0]
    - [1.0, 1.0, 1.0, 0.0]
  fitness_model:
    distribution: !Exponential
        weights:
          beneficial: 0.29
          deleterious: 0.51
          lethal: 0.2
          neutral: 0.0
        lambda_beneficial: 0.03
        lambda_deleterious: 0.21
    utility: Linear
simulation_plan:
  - generation: "(25 + {}) % 50"
    event: transmission
    value: migration_fwd
  - generation: "{} % 50"
    event: transmission
    value: migration_rev
  - generation: "{} % 200"
    event: sample
    value: 1000
```
