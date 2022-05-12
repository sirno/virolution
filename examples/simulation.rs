extern crate virolution;

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

    let fitness_table = FitnessTable::new(&sequence, &4, distribution);

    let wt = Wildtype::create_wildtype(sequence);
    let ht = wt.borrow_mut().create_descendant(2, 0x01);
    let ht2 = ht.borrow_mut().create_descendant(1, 0x02);
    let ht3 = ht2.borrow_mut().create_descendant(2, 0x03);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2);

    println!("---fitnesses---");
    println!("wt: {}", wt.borrow().get_fitness(&fitness_table));
    println!("ht: {}", ht.borrow().get_fitness(&fitness_table));
    println!("ht2: {}", ht2.borrow().get_fitness(&fitness_table));
    println!("ht3: {}", ht3.borrow().get_fitness(&fitness_table));
    println!("ht4: {}", ht4.borrow().get_fitness(&fitness_table));

    let population = (0..10).map(|_| wt.get_clone()).collect();
    let simulation_settings = SimulationSettings {
        mutation_rate: 1e-2,
        substitution_matrix: [
            [0., 1., 1., 1.],
            [1., 0., 1., 1.],
            [1., 1., 0., 1.],
            [1., 1., 1., 0.],
        ],
        host_population_size: 5,
        infection_fraction: 0.7,
        basic_reproductive_number: 100.,
        max_population: 100,
        dilution: 0.17,
    };
    simulation_settings.write("settings.yaml");
    let settings = SimulationSettings::read("settings.yaml");
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
            .map(|hap| hap.borrow().get_string())
            .collect::<Vec<String>>()
    );
    println!("offspring: {:?}", offspring);
    println!(
        "subsample: {:?}",
        population2
            .iter()
            .map(|hap| hap.borrow().get_string())
            .collect::<Vec<String>>()
    );
}
