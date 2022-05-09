#![feature(let_chains)]
#![feature(test)]

extern crate test;

mod fitness;
mod haplotype;
mod simulation;
mod simulation_settings;
mod transfers;

use fitness::*;
use haplotype::*;
use simulation::*;
use simulation_settings::*;
use std::rc::Rc;
use transfers::*;

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
    let init_population = (0..1_000_000).map(|_| Rc::clone(&wt)).collect();
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

fn _simulation_compartment_experiments() {
    let plan = Plan::read("plan.csv");
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
    let settings = SimulationSettings {
        mutation_rate: 1e-6,
        substitution_matrix: [
            [0., 1., 1., 1.],
            [1., 0., 1., 1.],
            [1., 1., 0., 1.],
            [1., 1., 1., 0.],
        ],
        host_population_size: 100_000_00,
        infection_fraction: 0.7,
        basic_reproductive_number: 100.,
        max_population: 100_000,
        dilution: 0.02,
    };

    let n_compartments = 3;
    let mut compartment_simulations: Vec<Simulation> = (0..n_compartments)
        .map(|_| {
            let init_population = (0..100_000).map(|_| Rc::clone(&wt)).collect();
            Simulation::new(
                Rc::clone(&wt),
                init_population,
                fitness_table.clone(),
                settings.clone(),
            )
        })
        .collect();

    for gen in 1..=20 {
        println!("generation={}", gen);
        let mut offsprings: Vec<Vec<usize>> = Vec::new();
        for simulation in compartment_simulations.iter_mut() {
            let infectant_map = simulation.get_infectant_map();
            let host_map = simulation.get_host_map(&infectant_map);
            simulation.mutate_infectants(&host_map);
            let offspring = simulation.replicate_infectants(&host_map);
            offsprings.push(offspring);
        }
        let mut populations: Vec<Population> = vec![Vec::new(); n_compartments];
        let transfers = plan.get_transfer_matrix(gen);
        for origin in 0..n_compartments {
            for target in 0..n_compartments {
                let mut population = compartment_simulations[origin]
                    .subsample_population(&offsprings[origin], transfers[origin][target]);
                populations[target].append(&mut population);
            }
        }
        for (idx, simulation) in compartment_simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone());
        }
    }
}

fn main() {
    _haplotype_experiments();
    _population_experiments();
    _simulation_experiments();
    // _simulation_loop_experiments();
    _simulation_compartment_experiments();
}

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
        let init_population: Population = (0..10).map(|_| Rc::clone(&wt)).collect();
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
