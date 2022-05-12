#![feature(let_chains)]
#![feature(test)]

extern crate test;

use rayon::prelude::*;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::simulation::*;
use virolution::simulation_settings::*;
use virolution::transfers::*;

fn _simulation_loop_experiments() {
    let sequence = vec![Some(0x00); 5386];
    let distribution = FitnessDistribution::Exponential(ExponentialParameters {
        weights: MutationCategoryWeights {
            beneficial: 0.29,
            deleterious: 0.51,
            lethal: 0.2,
            neutral: 0.,
        },
        lambda_beneficial: 0.03,
        lambda_deleterious: 0.21,
    });

    let fitness_table = FitnessTable::new(&sequence, &4, distribution);

    let wt = Wildtype::create_wildtype(sequence);
    let init_population = (0..1_000_000).map(|_| wt.get_clone()).collect();
    let settings = SimulationSettings {
        mutation_rate: 1e-6,
        substitution_matrix: [
            [0., 1., 1., 1.],
            [1., 0., 1., 1.],
            [1., 1., 0., 1.],
            [1., 1., 1., 0.],
        ],
        host_population_size: 1_000_000,
        infection_fraction: 0.7,
        basic_reproductive_number: 100.,
        max_population: 1_000_000,
        dilution: 0.02,
    };

    let mut simulation = Simulation::new(wt, init_population, fitness_table, settings);
    for gen in 1..=5 {
        println!("generation={}", gen);
        let infectant_map = simulation.get_infectant_map();
        let host_map = simulation.get_host_map(&infectant_map);
        simulation.mutate_infectants(&host_map);
        let offspring = simulation.replicate_infectants(&host_map);
        let population = simulation.subsample_population(&offspring, 1.);
        println!("pop_size: {}", population.len());
        simulation.set_population(population);
        // println!("{:?}", simulation.print_population());
    }
}

fn main() {}

#[cfg(test)]
mod tests {
    use super::*;
    use test::Bencher;

    #[bench]
    fn bench_next_generation(b: &mut Bencher) {
        let sequence = vec![Some(0x00); 100];
        let distribution = FitnessDistribution::Exponential(ExponentialParameters {
            weights: MutationCategoryWeights {
                beneficial: 0.29,
                deleterious: 0.51,
                lethal: 0.2,
                neutral: 0.,
            },
            lambda_beneficial: 0.03,
            lambda_deleterious: 0.21,
        });

        let fitness_table = FitnessTable::new(&sequence, &4, distribution);

        let wt = Wildtype::create_wildtype(sequence);
        let init_population: Population = (0..10).map(|_| wt.get_clone()).collect();
        let settings = SimulationSettings {
            mutation_rate: 1e-3,
            substitution_matrix: [
                [0., 1., 1., 1.],
                [1., 0., 1., 1.],
                [1., 1., 0., 1.],
                [1., 1., 1., 0.],
            ],
            host_population_size: 5,
            infection_fraction: 0.7,
            basic_reproductive_number: 100.,
            max_population: 100000000,
            dilution: 0.17,
        };
        let mut simulation = Simulation::new(wt, init_population, fitness_table, settings);
        b.iter(|| {
            simulation.next_generation();
        })
    }
}
