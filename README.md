# Virolation2000

This tool simulates a virus population under static selective pressure with
single-nucleotide polymorphism and recombination in discrete generations.
The host population is kept constant for each generation, and multiple
simulations can be run in parallel with controlled mixing. As input, the
simulation requires the model parameters, a wildtype sequence and a file with
all migration and sampling events. As output, a fasta file with samples of the
populations is generated.

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
pub struct SimulationSettings {
    pub basic_reproductive_number: f64,
    pub mutation_rate: f64,
    pub substitution_matrix: [[f64; 4]; 4],
    pub recombination_rate: f64,
    pub host_population_size: usize,
    pub max_population: usize,
    pub infection_fraction: f64,
    pub dilution: f64,
}
```
