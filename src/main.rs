#![feature(let_chains)]
#![feature(once_cell)]
#![feature(test)]

extern crate test;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use seq_io::fasta;
use seq_io::fasta::Record;
use std::fs;
use std::io;
use std::panic::catch_unwind;
use std::path::Path;
use virolution::args::*;
use virolution::barcode::*;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::references::HaplotypeRef;
use virolution::simulation::*;
use virolution::simulation_settings::*;
use virolution::transfers::*;

fn setup(args: &Args) {
    // setup logger
    simple_logging::log_to_file(args.log_file.as_str(), log::LevelFilter::Info)
        .expect("Failed to init logging.");

    // setup barcode file
    let barcode_path = Path::new(&args.outdir).join("barcodes.csv");
    std::fs::create_dir_all(barcode_path.parent().unwrap()).expect("Unable to create output path.");
    let mut barcode_file = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(barcode_path)
        .expect("Unable to open barcode file.");
    BarcodeEntry::write_header(&mut barcode_file).expect("Unable to write barcode header.");
}

fn load_sequence(path: &str) -> Vec<Symbol> {
    let mut reader = fasta::Reader::from_path(path).expect("Unable to open sequence file.");
    reader
        .next()
        .unwrap()
        .expect("Unable to read sequence.")
        .seq()
        .iter()
        .filter(|&&enc| enc != 0x0au8)
        .map(|enc| match catch_unwind(|| Some(FASTA_DECODE[enc])) {
            Ok(result) => result,
            Err(_) => panic!("Unable to decode literal {enc}."),
        })
        .collect()
}

#[cfg(feature = "parallel")]
fn create_simulations(
    args: &Args,
    wildtype: &HaplotypeRef,
    fitness_table: FitnessTable,
    settings: SimulationSettings,
) -> Vec<Simulation> {
    (0..args.n_compartments)
        .into_par_iter()
        .map(|compartment_idx| {
            let init_population: Population = if compartment_idx == 0 {
                (0..args.initial_population_size)
                    .into_par_iter()
                    .map(|_| wildtype.get_clone())
                    .collect()
            } else {
                Vec::new()
            };
            Simulation::new(
                wildtype.get_clone(),
                init_population,
                fitness_table.clone(),
                settings.clone(),
            )
        })
        .collect()
}

#[cfg(not(feature = "parallel"))]
fn create_simulations(
    args: &Args,
    wildtype: &HaplotypeRef,
    fitness_table: FitnessTable,
    settings: SimulationSettings,
) -> Vec<Simulation> {
    (0..args.n_compartments)
        .into_iter()
        .map(|compartment_idx| {
            let init_population: Population = if compartment_idx == 0 {
                (0..args.initial_population_size)
                    .into_iter()
                    .map(|_| wildtype.get_clone())
                    .collect()
            } else {
                Vec::new()
            };
            Simulation::new(
                wildtype.get_clone(),
                init_population,
                fitness_table.clone(),
                settings.clone(),
            )
        })
        .collect()
}

#[cfg(feature = "parallel")]
fn run(args: &Args, simulations: &mut Vec<Simulation>, plan: Plan) {
    // init progress bar
    let bar = ProgressBar::new(args.generations as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("[{bar:40}] {pos:>7}/{len:7} [{elapsed_precise} / {duration_precise}] {msg}")
            .expect("Unable to create template.")
            .progress_chars("=> "),
    );

    for generation in 0..=args.generations {
        // simulate compartmentalized population in parallel
        let offsprings: Vec<Vec<usize>> = simulations
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
            .collect();

        // transfer between compartments
        let transfer = plan.get_transfer_matrix(generation);

        let populations: Vec<Population> = (0..args.n_compartments)
            .into_par_iter()
            .map(|target| {
                let target_populations: Vec<Population> = (0..args.n_compartments)
                    .into_iter()
                    .map(|origin| {
                        simulations[origin]
                            .subsample_population(&offsprings[origin], transfer[target][origin])
                    })
                    .collect();

                target_populations.concat()
            })
            .collect();

        // update populations
        for (idx, simulation) in simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone());
        }

        // logging
        let population_sizes: Vec<usize> = populations.iter().map(|pop| pop.len()).collect();

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );

        // write to output when sampling
        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            log::info!("Sampling {} individuals...", sample_size);
            for (compartment_id, compartment) in simulations.iter().enumerate() {
                let barcode = format!("sample_{generation}_{compartment_id}");

                // create output files
                let barcode_path = Path::new(&args.outdir).join("barcodes.csv");
                let sample_path = Path::new(&args.outdir).join(format!("{barcode}.fasta"));

                // create file buffers
                let mut barcode_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(barcode_path)
                    .expect("Unable to open barcode file.");
                let mut samples_file = io::BufWriter::new(fs::File::create(sample_path).unwrap());

                // sample sequences and write to file
                for (sequence_id, sequence) in compartment
                    .get_population()
                    .choose_multiple(&mut rand::thread_rng(), sample_size)
                    .enumerate()
                {
                    let record = sequence.get_record(
                        format!(
                            "compartment_id={};sequence_id={};generation={}",
                            compartment_id, sequence_id, generation
                        )
                        .as_str(),
                    );
                    record
                        .write(&mut samples_file)
                        .expect("Unable to write to file.");
                }

                // write barcode to file
                BarcodeEntry {
                    barcode: &barcode,
                    experiment: &args.simulation_name,
                    time: generation,
                    replicate: 0,
                    compartment: compartment_id,
                }
                .write(&mut barcode_file)
                .expect("Unable to write to barcode file.");
            }
        }

        // update progress bar
        bar.set_position(generation.try_into().unwrap());
        bar.set_message(format!("{population_sizes:?}"));
    }
    bar.finish_with_message("Done.");
    log::info!("Finished simulation.");
}

#[cfg(not(feature = "parallel"))]
fn run(args: &Args, simulations: &mut Vec<Simulation>, plan: Plan) {
    // init progress bar
    let bar = ProgressBar::new(args.generations as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("[{bar:40}] {pos:>7}/{len:7} [{elapsed_precise} / {duration_precise}] {msg}")
            .expect("Unable to create template.")
            .progress_chars("=> "),
    );

    for generation in 0..=args.generations {
        // simulate compartmentalized population in parallel
        let offsprings: Vec<Vec<usize>> = simulations
            .iter_mut()
            .map(|simulation| {
                if simulation.get_population().is_empty() {
                    return Vec::new();
                }
                let infectant_map = simulation.get_infectant_map();
                let host_map = simulation.get_host_map(&infectant_map);
                simulation.mutate_infectants(&host_map);
                simulation.replicate_infectants(&host_map)
            })
            .collect();

        // transfer between compartments
        let transfer = plan.get_transfer_matrix(generation);

        let populations: Vec<Population> = (0..args.n_compartments)
            .into_iter()
            .map(|target| {
                let target_populations: Vec<Population> = (0..args.n_compartments)
                    .into_iter()
                    .map(|origin| {
                        simulations[origin]
                            .subsample_population(&offsprings[origin], transfer[target][origin])
                    })
                    .collect();

                target_populations.concat()
            })
            .collect();

        // update populations
        for (idx, simulation) in simulations.iter_mut().enumerate() {
            simulation.set_population(populations[idx].clone());
        }

        // logging
        let population_sizes: Vec<usize> = populations.iter().map(|pop| pop.len()).collect();

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );

        // write to output when sampling
        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            log::info!("Sampling {} individuals...", sample_size);
            for (compartment_id, compartment) in simulations.iter().enumerate() {
                let barcode = format!("sample_{generation}_{compartment_id}");

                // create output files
                let barcode_path = Path::new(&args.outdir).join("barcodes.csv");
                let sample_path = Path::new(&args.outdir).join(format!("{barcode}.fasta"));

                // create file buffers
                let mut barcode_file = fs::OpenOptions::new()
                    .create(true)
                    .append(true)
                    .open(barcode_path)
                    .expect("Unable to open barcode file.");
                let mut samples_file = io::BufWriter::new(fs::File::create(sample_path).unwrap());

                // sample sequences and write to file
                for (sequence_id, sequence) in compartment
                    .get_population()
                    .choose_multiple(&mut rand::thread_rng(), sample_size)
                    .enumerate()
                {
                    let record = sequence.get_record(
                        format!(
                            "compartment_id={};sequence_id={};generation={}",
                            compartment_id, sequence_id, generation
                        )
                        .as_str(),
                    );
                    record
                        .write(&mut samples_file)
                        .expect("Unable to write to file.");
                }

                // write barcode to file
                BarcodeEntry {
                    barcode: &barcode,
                    experiment: &args.simulation_name,
                    time: generation,
                    replicate: 0,
                    compartment: compartment_id,
                }
                .write(&mut barcode_file)
                .expect("Unable to write to barcode file.");
            }
        }

        // update progress bar
        bar.set_position(generation.try_into().unwrap());
        bar.set_message(format!("{population_sizes:?}"));
    }
    bar.finish_with_message("Done.");
    log::info!("Finished simulation.");
}

fn main() {
    if cfg!(feature = "parallel") {
        println!("Running in parallel mode.");
    } else {
        println!("Running in serial mode.");
    }

    let args = Args::parse();

    setup(&args);

    let plan = Plan::read(args.transfer_plan.as_str());
    let sequence = load_sequence(args.sequence.as_str());

    // create wildtype
    let wildtype = Wildtype::new(sequence.clone());

    // load simulation settings
    let settings = SimulationSettings::read_from_file(args.settings.as_str()).unwrap();
    log::info!("Loaded settings\n{}", settings);

    // create and write fitness table
    let fitness_table = FitnessTable::new(&sequence, 4, settings.fitness_distribution.clone());
    let mut fitness_file =
        io::BufWriter::new(fs::File::create(args.fitness_table.clone()).unwrap());
    fitness_table.write(&mut fitness_file).unwrap();

    // create individual compartments
    println!("Creating {} compartments...", args.n_compartments);
    let mut simulations: Vec<Simulation> =
        create_simulations(&args, &wildtype, fitness_table, settings);

    // run simulation
    run(&args, &mut simulations, plan);

    // store tree if specified.
    log::info!("Storing tree...");
    if let Some(tree_file) = args.trees {
        fs::write(tree_file, wildtype.get_tree())
            .unwrap_or_else(|_| eprintln!("Unable to write tree file."));
    }
    log::info!("Finished storing tree.");
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

        let fitness_table = FitnessTable::new(&sequence, 4, distribution.clone());

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
            fitness_distribution: distribution,
        };
        let mut simulation = Simulation::new(wt, init_population, fitness_table, settings);
        b.iter(|| {
            simulation.next_generation();
        })
    }
}
