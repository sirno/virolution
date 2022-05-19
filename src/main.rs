#![feature(let_chains)]
#![feature(test)]

extern crate test;

use clap::Parser;
use rand::prelude::*;
use rayon::prelude::*;
use seq_io::fasta;
use seq_io::fasta::Record;
use std::fs;
use std::io;
use virolution::args::*;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::simulation::*;
use virolution::simulation_settings::*;
use virolution::transfers::*;

fn main() {
    let args = Args::parse();
    let plan = Plan::read(args.transfer_plan.as_str());
    let mut reader =
        fasta::Reader::from_path(args.sequence).expect("Unable to open sequence file.");
    let sequence: Vec<Option<u8>> = reader
        .next()
        .unwrap()
        .expect("Unable to read sequence.")
        .seq()
        .iter()
        .map(|enc| Some(FASTA_DECODE[enc]))
        .collect();
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
    let settings = SimulationSettings::read(args.settings.as_str());

    let n_compartments = 3;
    let mut compartment_simulations: Vec<Simulation> = (0..n_compartments)
        .map(|_| {
            let init_population = (0..1_000_000).map(|_| wt.get_clone()).collect();
            Simulation::new(
                wt.get_clone(),
                init_population,
                fitness_table.clone(),
                settings.clone(),
            )
        })
        .collect();

    let mut writer = io::BufWriter::new(fs::File::create(args.output).unwrap());

    for generation in 1..=args.generations {
        println!("generation={}", generation);

        let mut offsprings: Vec<Vec<usize>> = Vec::new();
        compartment_simulations
            .par_iter_mut()
            .map(|simulation| {
                let infectant_map = simulation.get_infectant_map();
                let host_map = simulation.get_host_map(&infectant_map);
                simulation.mutate_infectants(&host_map);
                simulation.replicate_infectants(&host_map)
            })
            .collect_into_vec(&mut offsprings);

        let mut populations: Vec<Population> = vec![Vec::new(); n_compartments];
        let transfers = plan.get_transfer_matrix(generation);
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

        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            for (compartment_id, compartment) in compartment_simulations.iter().enumerate() {
                for (sequence_id, sequence) in compartment
                    .get_population()
                    .choose_multiple(&mut rand::thread_rng(), sample_size)
                    .enumerate()
                {
                    let record = sequence.borrow().get_record(
                        format!(
                            "compartment_id={};sequence_id={};generation={}",
                            compartment_id, sequence_id, generation
                        )
                        .as_str(),
                    );
                    record.write(&mut writer).expect("Unable to write to file.");
                    if sequence_id > 10 {
                        break;
                    }
                }
            }
        }
    }
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
        let init_population: Population = (0..10).map(|_| wt.get_clone()).collect();
        let settings = SimulationSettings {
            mutation_rate: 1e-3,
            recombination_rate: 1e-5,
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
