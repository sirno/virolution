#![feature(let_chains)]
#![feature(test)]

extern crate test;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
use rayon::prelude::*;
use seq_io::fasta;
use seq_io::fasta::Record;
use std::fs;
use std::io;
use std::panic::catch_unwind;
use virolution::args::*;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::simulation::*;
use virolution::simulation_settings::*;
use virolution::transfers::*;

fn main() {
    let args = Args::parse();
    simple_logging::log_to_file(args.log_file.as_str(), log::LevelFilter::Info)
        .expect("Failed to init logging.");
    let plan = Plan::read(args.transfer_plan.as_str());
    let mut reader =
        fasta::Reader::from_path(args.sequence).expect("Unable to open sequence file.");
    let sequence: Vec<Option<u8>> = reader
        .next()
        .unwrap()
        .expect("Unable to read sequence.")
        .seq()
        .iter()
        .filter(|&&enc| enc != 0x0au8)
        .map(|enc| match catch_unwind(|| Some(FASTA_DECODE[enc])) {
            Ok(result) => result,
            Err(_) => panic!("Unable to decode literal {}.", enc),
        })
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

    let wt = Wildtype::new(sequence);
    let settings = SimulationSettings::read(args.settings.as_str());

    // create individual compartments
    let n_compartments = args.n_compartments;
    print!("Creating {} compartments...\n", n_compartments);
    let mut compartment_simulations: Vec<Simulation> = (0..n_compartments)
        .map(|compartment_idx| {
            let init_population: Population = if compartment_idx == 0 {
                (0..1_000_000).map(|_| wt.get_clone()).collect()
            } else {
                Vec::new()
            };
            Simulation::new(
                wt.get_clone(),
                init_population,
                fitness_table.clone(),
                settings.clone(),
            )
        })
        .collect();

    // create output files
    let mut writer = io::BufWriter::new(fs::File::create(args.output).unwrap());

    // init progress bar
    let bar = ProgressBar::new(args.generations as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("[{bar:40}] {pos:>7}/{len:7} [{elapsed_precise} / {duration_precise}] {msg}")
            .progress_chars("=> "),
    );

    // run simulation
    for generation in 0..=args.generations {
        // simulate compartmentalized population in parallel
        let mut offsprings: Vec<Vec<usize>> = Vec::new();
        compartment_simulations
            .par_iter_mut()
            .map(|simulation| {
                if simulation.get_population().is_empty() {
                    return Vec::new();
                }
                let infectant_map = simulation.get_infectant_map();
                let host_map = simulation.get_host_map(&infectant_map);
                simulation.mutate_infectants(&host_map);
                simulation.replicate_infectants(&host_map)
            })
            .collect_into_vec(&mut offsprings);

        // transfer between compartments
        let mut populations: Vec<Population> = vec![Vec::new(); n_compartments];
        let transfers = plan.get_transfer_matrix(generation);
        for origin in 0..n_compartments {
            #[allow(clippy::needless_range_loop)]
            for target in 0..n_compartments {
                let mut population = compartment_simulations[origin]
                    .subsample_population(&offsprings[origin], transfers[target][origin]);
                populations[target].append(&mut population);
            }
        }

        // write to output when sampling
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

        // update populations
        for (idx, simulation) in compartment_simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone());
        }

        // logging
        let population_sizes: Vec<usize> = populations.iter().map(|pop| pop.len()).collect();

        bar.set_position(generation.try_into().unwrap());
        bar.set_message(format!("{population_sizes:?}"));

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );
    }
    bar.finish_with_message("Done.");

    // Store tree if specified.
    if let Some(tree_file) = args.trees {
        fs::write(tree_file, wt.borrow().get_tree())
            .unwrap_or_else(|_| eprintln!("Unable to write tree file."));
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

        let wt = Wildtype::new(sequence);
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
