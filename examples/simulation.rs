extern crate virolution;

use std::fs;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::simulation::*;
use virolution::simulation_settings::*;

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

    let fitness_table = FitnessTable::new(&sequence, 4, distribution.clone());

    let wt = Wildtype::new(sequence);
    let ht = wt.create_descendant(2, 0x01);
    let ht2 = ht.create_descendant(1, 0x02);
    let ht3 = ht2.create_descendant(2, 0x03);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2);

    println!("---fitnesses---");
    println!("wt: {}", wt.get_fitness(&fitness_table));
    println!("ht: {}", ht.get_fitness(&fitness_table));
    println!("ht2: {}", ht2.get_fitness(&fitness_table));
    println!("ht3: {}", ht3.get_fitness(&fitness_table));
    println!("ht4: {}", ht4.get_fitness(&fitness_table));

    let population = (0..10).map(|_| wt.get_clone()).collect();
    let simulation_settings = SimulationSettings {
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
        fitness_distribution: distribution,
    };
    simulation_settings
        .write_to_file("settings_example.yaml")
        .expect("Failed to write settings to file");

    let settings = SimulationSettings::read_from_file("settings_example.yaml")
        .expect("Failed to read settings from file");
    fs::remove_file("settings_example.yaml").expect("Unable to remove file.");
    let mut simulation = Simulation::new(wt, population, fitness_table, settings.clone());
    let infectant_map = simulation.get_infectant_map();
    let host_map = simulation.get_host_map(&infectant_map);
    simulation.mutate_infectants(&host_map);
    let offspring = simulation.replicate_infectants(&host_map);
    let population2 = simulation.subsample_population(&offspring, 1.);

    println!("---simulation_settings---");
    println!("{:#?}", settings);
    println!("---infection-mapping---");
    println!("infectant_map: {:?}", infectant_map);
    println!("host_map: {:?}", host_map);
    println!(
        "population: {:?}",
        simulation
            .get_population()
            .iter()
            .map(|hap| hap.get_string())
            .collect::<Vec<String>>()
    );
    println!("offspring: {:?}", offspring);
    println!(
        "subsample: {:?}",
        population2
            .iter()
            .map(|hap| hap.get_string())
            .collect::<Vec<String>>()
    );
}
