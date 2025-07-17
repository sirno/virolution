extern crate virolution;

use std::path::PathBuf;
use std::sync::Arc;

use virolution::config::{HostFitness, HostModel, Parameters, Schedule};
use virolution::core::fitness::UtilityFunction;
use virolution::core::haplotype::*;
use virolution::core::population::{Population, Store};
use virolution::encoding::Nucleotide as Nt;
use virolution::init::fitness::*;
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
    let host_fitness = HostFitness::new(
        Some(FitnessModel::new(
            distribution.clone(),
            UtilityFunction::Linear,
        )),
        None,
    );

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
        basic_reproductive_number: 100.,
        max_population: 1_000_000,
        dilution: 0.02,
        host_model: HostModel::SingleHost(host_fitness),
        n_hits: 1,
    };

    let (mut attribute_definition, host_specs) =
        parameters
            .host_model
            .make_definitions(&parameters, &sequence, None);

    let generation_provider = Arc::new(Generation::new(0));
    attribute_definition.register(generation_provider.clone());

    let wt = Wildtype::new(sequence, &attribute_definition);

    let n_compartments = 3;
    let mut compartment_simulations: Vec<BasicSimulation<Nt>> = (0..n_compartments)
        .map(|_| {
            let population = population![wt.clone(), 1_000_000];
            BasicSimulation::new(
                wt.clone(),
                population,
                parameters.clone(),
                host_specs.clone(),
                generation_provider.clone(),
            )
        })
        .collect();

    for generation in 1..=50 {
        println!("generation={generation}");
        let offsprings: Vec<Vec<usize>> = compartment_simulations
            .iter_mut()
            .map(|simulation| {
                simulation.infect();
                simulation.mutate_infectants();
                simulation.replicate_infectants()
            })
            .collect();

        let transfers = plan.get_transfer_matrix(generation);

        let populations: Vec<Population<Store<Nt>>> = (0..n_compartments)
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
