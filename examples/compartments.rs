extern crate virolution;

use std::path::PathBuf;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::population;
use virolution::population::Population;
use virolution::simulation::*;
use virolution::simulation_parameters::*;
use virolution::simulation_plan::*;

fn main() {
    let plan_path = PathBuf::from_iter([env!("CARGO_MANIFEST_DIR"), "data/plan.csv"]);
    let plan = SimulationPlan::read(plan_path.to_str().unwrap()).expect("Failed to read plan");
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
    let fitness_model = FitnessModel::new(distribution.clone(), UtilityFunction::Linear);

    let fitness_table = FitnessTable::from_model(&sequence, 4, fitness_model.clone()).unwrap();

    let wt = Wildtype::new(sequence);
    let settings = SimulationParameters {
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
    let mut compartment_simulations: Vec<BasicSimulation> = (0..n_compartments)
        .map(|_| {
            let population = population![wt.clone(); 1_000_000];
            let fitness_tables = vec![(0..settings.host_population_size, fitness_table.clone())];
            BasicSimulation::new(
                wt.get_clone(),
                population,
                fitness_tables,
                settings.clone(),
                0,
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
                        .subsample_population(&offsprings[origin], *transfers.get(target, origin))
                }))
            })
            .collect();

        for (idx, simulation) in compartment_simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone());
        }
    }
}
