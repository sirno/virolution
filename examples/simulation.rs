extern crate virolution;

use std::fs;
use std::sync::Arc;
use virolution::config::{FitnessModelField, HostFitness, Parameters};
use virolution::core::attributes::AttributeSetDefinition;
use virolution::core::fitness::UtilityFunction;
use virolution::core::haplotype::*;
use virolution::core::hosts::{HostSpec, HostSpecs};
use virolution::core::population::Store;
use virolution::core::Population;
use virolution::encoding::Nucleotide as Nt;
use virolution::init::fitness::*;
use virolution::providers::FitnessProvider;
use virolution::providers::Generation;
use virolution::simulation::*;

use virolution::population;

fn main() {
    let sequence = vec![Nt::A; 100];
    let sequence_length = sequence.len();
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

    let mut attribute_definition = AttributeSetDefinition::new();
    let name = "fitness";
    attribute_definition.register(Arc::new(
        FitnessProvider::from_model(name, &sequence, &fitness_model).unwrap(),
    ));
    let generation_provider = Arc::new(Generation::new(0));
    attribute_definition.register(generation_provider.clone());
    let wt = Wildtype::new(sequence, &attribute_definition);
    let ht = wt.create_descendant(vec![2], vec![Nt::T]);
    let ht2 = ht.create_descendant(vec![1], vec![Nt::C]);
    let ht3 = ht2.create_descendant(vec![2], vec![Nt::G]);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2);

    println!("---fitnesses---");
    println!("wt: {}", wt.get_attribute_or_compute("fitness").unwrap());
    println!("ht: {}", ht.get_attribute_or_compute("fitness").unwrap());
    println!("ht2: {}", ht2.get_attribute_or_compute("fitness").unwrap());
    println!("ht3: {}", ht3.get_attribute_or_compute("fitness").unwrap());
    println!("ht4: {}", ht4.get_attribute_or_compute("fitness").unwrap());

    let population: Population<Store<Nt>> = population![wt.clone(), 10];
    let simulation_settings = Parameters {
        mutation_rate: 1e-6,
        recombination_rate: 0.,
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
        fitness_model: FitnessModelField::SingleHost(HostFitness::new(Some(fitness_model), None)),
    };
    simulation_settings
        .write_to_file("parameters_example.yaml")
        .expect("Failed to write settings to file");

    let settings = Parameters::read_from_file("parameters_example.yaml")
        .expect("Failed to read settings from file");
    fs::remove_file("parameters_example.yaml").expect("Unable to remove file.");
    let host_specs = HostSpecs::from_vec(vec![HostSpec::new(
        0..simulation_settings.host_population_size,
        Box::new(BasicHost::new(
            sequence_length,
            &simulation_settings,
            Some(name),
            None,
        )),
    )]);
    let mut simulation = BasicSimulation::new(
        wt,
        population,
        settings.clone(),
        host_specs,
        generation_provider,
    );
    simulation.infect();
    simulation.mutate_infectants();
    let offspring = simulation.replicate_infectants();
    let population2 = simulation.subsample_population(&offspring, 1.);

    println!("---simulation_settings---");
    println!("{:#?}", settings);
    println!(
        "population: {:?}",
        simulation
            .get_population()
            .iter()
            .map(|haplotype| haplotype.get_string())
            .collect::<Vec<String>>()
    );
    println!("offspring: {:?}", offspring);
    println!(
        "subsample: {:?}",
        population2
            .iter()
            .map(|haplotype| haplotype.get_string())
            .collect::<Vec<String>>()
    );
}
