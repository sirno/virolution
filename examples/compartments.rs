extern crate virolution;

use std::collections::HashMap;
use std::path::PathBuf;
use virolution::fitness::*;
use virolution::haplotype::*;
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

    let fitness_table = FitnessTable::new(&sequence, 4, fitness_model);

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
        fitness_distribution: distribution,
    };

    let n_compartments = 3;
    let mut compartment_simulations: Vec<Simulation> = (0..n_compartments)
        .map(|_| {
            let genotypes = HashMap::from_iter([(wt.get_id(), wt.get_clone())]);
            let init_population = (0..1_000_000).map(|_| wt.get_id()).collect();
            Simulation::new(
                wt.get_clone(),
                init_population,
                genotypes,
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
                let infectant_map = simulation.get_infectant_map();
                let host_map = simulation.get_host_map(&infectant_map);
                simulation.mutate_infectants(&host_map);
                simulation.replicate_infectants(&host_map)
            })
            .collect();

        let mut populations: Vec<Population> = vec![Vec::new(); n_compartments];
        let mut genotypes: Vec<Genotypes> = vec![HashMap::new(); n_compartments];
        let transfers = plan.get_transfer_matrix(gen);
        for origin in 0..n_compartments {
            for target in 0..n_compartments {
                if transfers[origin][target] == 0. {
                    continue;
                }
                let mut population = compartment_simulations[origin]
                    .subsample_population(&offsprings[origin], transfers[origin][target]);
                populations[target].append(&mut population);
                genotypes[target].extend(compartment_simulations[origin].get_genotypes());
            }
        }

        for (idx, simulation) in compartment_simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone(), Some(&genotypes[idx]));
        }
    }
}
