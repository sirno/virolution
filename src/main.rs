#![feature(let_chains)]

mod fitness;
mod haplotype;
mod simulation;

use fitness::*;
use haplotype::*;
use simulation::*;
use std::rc::Rc;

fn _haplotype_experiments() {
    let bytes = vec![Some(0x00); 4];
    let wt = Wildtype::create_wildtype(bytes);
    let ht = wt.borrow_mut().create_descendant(2, 0x01);
    let ht2 = ht.borrow_mut().create_descendant(1, 0x02);
    let ht3 = ht2.borrow_mut().create_descendant(2, 0x03);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2);
    println!("---debug---");
    println!("wt: {:?}", *wt.borrow());
    println!("ht: {:?}", *ht.borrow());
    println!("ht2: {:?}", *ht2.borrow());
    println!("ht3: {:?}", *ht3.borrow());
    println!("ht4: {:?}", *ht4.borrow());
    println!("---get_base---");
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| wt.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht2.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht3.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht4.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );

    println!("---get_sequence---");
    println!("{:?}", wt.borrow().get_sequence());
    println!("{:?}", ht.borrow().get_sequence());
    println!("{:?}", ht2.borrow().get_sequence());
    println!("{:?}", ht3.borrow().get_sequence());
    println!("{:?}", ht4.borrow().get_sequence());
}

fn _population_experiments() {}

fn _simulation_experiments() {
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

    let population = (0..10).map(|_| Rc::clone(&wt)).collect();
    let mut simulation = Simulation::new(population, fitness_table, 5, 0.7, 100., 1e-2);
    let infectant_map = simulation.get_infectant_map();
    let host_map = simulation.get_host_map(&infectant_map);
    simulation.mutate_infectants(&host_map);
    let offspring = simulation.replicate_infectants(&host_map);
    let population2 = simulation.subsample_population(&offspring);

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
    )
}

fn main() {
    _haplotype_experiments();
    _population_experiments();
    _simulation_experiments();
}
