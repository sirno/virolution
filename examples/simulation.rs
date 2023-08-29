extern crate virolution;

use std::fs;
use virolution::config::{FitnessModelField, Parameters};
use virolution::core::fitness::*;
use virolution::core::haplotype::*;
use virolution::core::Population;
use virolution::simulation::*;

use virolution::population;

fn main() {
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
    let fitness_model = FitnessModel::new(distribution.clone(), UtilityFunction::Linear);

    let fitness_table = FitnessTable::from_model(&sequence, 4, fitness_model.clone()).unwrap();

    let wt = Wildtype::new(sequence);
    let ht = wt.create_descendant(vec![2], vec![Some(0x01)], 0);
    let ht2 = ht.create_descendant(vec![1], vec![Some(0x02)], 0);
    let ht3 = ht2.create_descendant(vec![2], vec![Some(0x03)], 0);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2, 0);

    println!("---fitnesses---");
    println!("wt: {}", wt.get_fitness(&fitness_table));
    println!("ht: {}", ht.get_fitness(&fitness_table));
    println!("ht2: {}", ht2.get_fitness(&fitness_table));
    println!("ht3: {}", ht3.get_fitness(&fitness_table));
    println!("ht4: {}", ht4.get_fitness(&fitness_table));

    let population: Population = population![wt.clone(); 10];
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
        fitness_model: FitnessModelField::SingleHost(fitness_model),
    };
    simulation_settings
        .write_to_file("parameters_example.yaml")
        .expect("Failed to write settings to file");

    let settings = Parameters::read_from_file("parameters_example.yaml")
        .expect("Failed to read settings from file");
    fs::remove_file("parameters_example.yaml").expect("Unable to remove file.");
    let fitness_tables = vec![(0..simulation_settings.host_population_size, fitness_table)];
    let mut simulation = BasicSimulation::new(wt, population, fitness_tables, settings.clone(), 0);
    let host_map = simulation.get_host_map();
    simulation.mutate_infectants(&host_map);
    let offspring = simulation.replicate_infectants(&host_map);
    let population2 = simulation.subsample_population(&offspring, 1.);

    println!("---simulation_settings---");
    println!("{:#?}", settings);
    println!("---infection-mapping---");
    println!("host_map: {:?}", host_map);
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
