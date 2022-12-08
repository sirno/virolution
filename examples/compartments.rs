extern crate virolution;

use std::path::PathBuf;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::population;
use virolution::population::Population;
use virolution::simulation::*;
use virolution::simulation_settings::*;
use virolution::transfers::*;

fn main() {
    let plan_path = PathBuf::from_iter([env!("CARGO_MANIFEST_DIR"), "data/plan.csv"]);
    let plan = Plan::read(plan_path.to_str().unwrap());
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
    let fitness_model = FitnessModel {
        distribution: distribution.clone(),
        utility: UtilityFunction::Linear,
    };

    let fitness_table = FitnessTable::new(&sequence, 4, fitness_model.clone());

    let wt = Wildtype::new(sequence);
    let settings = SimulationSettings {
        mutation_rate: 1e-6,
        recombination_rate: 1e-8,
        substitution_matrix: [
            [0., 1., 1., 1.],
            [1., 0., 1., 1.],
            [1., 1., 0., 1.],
            [1., 1., 1., 0.],
        ],
        host_population_size: 10_000_000,
        infection_fraction: 0.7,
        basic_reproductive_number: 100.,
        max_population: 1_000_000,
        dilution: 0.02,
        fitness_model,
    };

    let n_compartments = 3;
    let mut compartment_simulations: Vec<Simulation> = (0..n_compartments)
        .map(|_| {
            let population = population![wt.clone(); 1_000_000];
            Simulation::new(
                wt.get_clone(),
                population,
                fitness_table.clone(),
                settings.clone(),
            )
        })
        .collect();

    for gen in 1..=50 {
        println!("generation={}", gen);
        let offsprings: Vec<Vec<usize>> = compartment_simulations
            .iter_mut()
            .map(|simulation| {
                let host_map = simulation.get_host_map();
                simulation.mutate_infectants(&host_map);
                simulation.replicate_infectants(&host_map)
            })
            .collect();

        let transfers = plan.get_transfer_matrix(gen);

        let populations: Vec<Population> = (0..n_compartments)
            .into_iter()
            .map(|target| {
                Population::from_iter((0..n_compartments).map(|origin| {
                    compartment_simulations[origin]
                        .subsample_population(&offsprings[origin], transfers[target][origin])
                }))
            })
            .collect();

        for (idx, simulation) in compartment_simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone());
        }
    }
}
