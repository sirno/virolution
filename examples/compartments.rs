extern crate virolution;

use std::path::PathBuf;
use std::sync::Arc;

use virolution::config::{FitnessModelField, Parameters, Schedule};
use virolution::core::attributes::AttributeSetDefinition;
use virolution::core::fitness::init::*;
use virolution::core::fitness::utility::UtilityFunction;
use virolution::core::fitness::FitnessProvider;
use virolution::core::haplotype::*;
use virolution::core::population::Population;
use virolution::encoding::Nucleotide as Nt;
use virolution::providers::Generation;
use virolution::simulation::*;

use virolution::population;

fn main() {
    let plan_path = PathBuf::from_iter([env!("CARGO_MANIFEST_DIR"), "data/plan.csv"]);
    let plan = Schedule::read(plan_path.to_str().unwrap()).expect("Failed to read plan");
    let sequence = vec![Nt::A; 5386];
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

    let name = "fitness";
    let fitness_provider = FitnessProvider::from_model(name, &sequence, &fitness_model)
        .expect("Failed to create fitness table");

    let generation_provider = Arc::new(Generation::new(0));

    let mut attribute_definition = AttributeSetDefinition::new();
    attribute_definition.register(Arc::new(fitness_provider));
    attribute_definition.register(generation_provider.clone());

    let wt = Wildtype::new(sequence, &attribute_definition);
    let parameters = Parameters {
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
        fitness_model: FitnessModelField::SingleHost(fitness_model),
    };

    let n_compartments = 3;
    let mut compartment_simulations: Vec<BasicSimulation<Nt>> = (0..n_compartments)
        .map(|_| {
            let population = population![wt.clone(), 1_000_000];
            let hosts = vec![(0..parameters.host_population_size, name)];
            BasicSimulation::new(
                wt.clone(),
                population,
                hosts,
                parameters.clone(),
                generation_provider.clone(),
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

        let populations: Vec<Population<Nt>> = (0..n_compartments)
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
