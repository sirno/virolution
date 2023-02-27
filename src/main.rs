#![feature(let_chains)]
#![feature(once_cell)]
#![feature(test)]

extern crate test;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use seq_io::fasta;
use seq_io::fasta::OwnedRecord;
use seq_io::fasta::Record;
use std::collections::HashMap;
use std::fs;
use std::io;
use std::panic::catch_unwind;
use std::path::Path;
use virolution::args::*;
use virolution::barcode::*;
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::plan::*;
use virolution::population;
use virolution::population::Population;
use virolution::references::HaplotypeRef;
use virolution::simulation::*;
use virolution::simulation_settings::*;

fn setup(args: &Args) {
    // setup logger
    let log_level = match args.verbose {
        0 => log::LevelFilter::Info,
        1 => log::LevelFilter::Debug,
        _ => log::LevelFilter::Trace,
    };
    simple_logging::log_to_file(args.log_file.as_str(), log_level).unwrap_or_else(|_| {
        eprintln!("Unable to open log file.");
        std::process::exit(1);
    });

    // setup rayon
    #[cfg(feature = "parallel")]
    if let Some(n_threads) = args.threads {
        println!("Setting number of threads to {}.", n_threads);
        rayon::ThreadPoolBuilder::new()
            .num_threads(n_threads)
            .build_global()
            .unwrap_or_else(|_| {
                eprintln!("Unable to set number of threads.");
                std::process::exit(1);
            });
    }

    // setup barcode file
    let barcode_path = Path::new(&args.outdir).join("barcodes.csv");
    std::fs::create_dir_all(barcode_path.parent().unwrap()).unwrap_or_else(|_| {
        eprintln!("Unable to create output path.");
        std::process::exit(1);
    });
    let mut barcode_file = fs::OpenOptions::new()
        .create(true)
        .append(true)
        .open(barcode_path)
        .expect("Unable to open barcode file.");
    BarcodeEntry::write_header(&mut barcode_file).expect("Unable to write barcode header.");
}

fn read_transfer_plan(path: &str) -> Plan {
    match Plan::read(path) {
        Ok(plan) => plan,
        Err(err) => match err {
            PlanReadError::IoError(err) => {
                eprintln!("Unable to open transfer plan from file '{path}'.");
                eprintln!("Reason: {err}");
                std::process::exit(1);
            }
            PlanReadError::CsvError(err) => {
                eprintln!("Unable to load transfer plan from file '{path}'.");
                eprintln!("Reason: {err}");
                std::process::exit(1);
            }
        },
    }
}

fn load_sequence(path: &str) -> Vec<Symbol> {
    let mut reader = fasta::Reader::from_path(path).unwrap_or_else(|err| match err.kind() {
        io::ErrorKind::NotFound => {
            eprintln!("Unable to find file {path}.");
            std::process::exit(1);
        }
        err => panic!("Unable to open file: {err}."),
    });
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

fn create_simulations(
    args: &Args,
    wildtype: &HaplotypeRef,
    fitness_table: FitnessTable,
    fitness_table2: FitnessTable,
    settings: SimulationSettings,
) -> Vec<BasicSimulation> {
    (0..args.n_compartments)
        .map(|compartment_idx| {
            let init_population: Population = if compartment_idx == 0 {
                population![wildtype.clone(); args.initial_population_size]
            } else {
                Population::new()
            };
            let fitness_tables = vec![
                (
                    0..(settings.host_population_size / 2),
                    fitness_table.clone(),
                ),
                (
                    ((settings.host_population_size / 2)..settings.host_population_size),
                    fitness_table2.clone(),
                ),
            ];
            BasicSimulation::new(
                wildtype.get_clone(),
                init_population,
                fitness_tables,
                settings.clone(),
                0,
            )
        })
        .collect()
}

fn sample(simulations: &[BasicSimulation], sample_size: usize, generation: usize, args: &Args) {
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

        // precompute sequences
        let population = compartment.get_population();
        let sequences: HashMap<usize, Vec<u8>> = population
            .iter()
            .unique()
            .map(|haplotype_ref| {
                let sequence = haplotype_ref
                    .get_sequence()
                    .into_iter()
                    .map(|symbol| match symbol {
                        Some(s) => FASTA_ENCODE[&s],
                        None => 0x2d,
                    })
                    .collect();
                (haplotype_ref.get_id(), sequence)
            })
            .collect();

        // sample sequences and write to file
        for (haplotype_id, haplotype_ref) in population
            .choose_multiple(&mut rand::thread_rng(), sample_size)
            .into_iter()
            .enumerate()
        {
            let head = format!(
                "compartment_id={};sequence_id={};generation={}",
                compartment_id, haplotype_id, generation
            )
            .as_bytes()
            .to_vec();
            let sequence = sequences[&haplotype_ref.get_id()].clone();
            let record = OwnedRecord {
                head,
                seq: sequence,
            };
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

#[cfg(feature = "parallel")]
fn run(args: &Args, simulations: &mut Vec<BasicSimulation>, plan: Plan) {
    let bar = match args.disable_progress_bar {
        true => None,
        false => {
            let bar = ProgressBar::new(args.generations as u64);
            bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "[{bar:40}] {pos:>7}/{len:7} [{elapsed_precise} / {duration_precise}] {msg}",
                )
                .expect("Unable to create template.")
                .progress_chars("=> "),
        );
            Some(bar)
        }
    };

    for generation in 0..=args.generations {
        // logging
        log::debug!("Generate logging message for generation {generation}...");
        let population_sizes: Vec<usize> = simulations
            .iter()
            .map(|sim| sim.get_population().len())
            .collect();

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );

        if let Some(bar) = bar.as_ref() {
            bar.set_position(generation.try_into().unwrap());
            bar.set_message(format!("{population_sizes:?}"));
        }

        // write to output when sampling
        log::debug!("Process sampling...");
        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            sample(simulations, sample_size, generation, args);
        }

        // abort on last generation after sampling
        if generation == args.generations {
            break;
        }

        // adjust settings if needed
        if let Some(settings) = plan.get_settings(generation) {
            log::info!("Adjusting settings to:\n{}", settings);
            simulations.iter_mut().for_each(|simulation| {
                simulation.set_settings(settings.clone());
            });
        }

        // increment generation
        simulations.iter_mut().for_each(|simulation| {
            simulation.increment_generation();
        });

        // simulate compartmentalized population in parallel
        log::debug!("Generate host maps...");
        let host_maps: Vec<HostMap> = simulations
            .par_iter()
            .map(|simulation| simulation.get_host_map())
            .collect();
        log::debug!("Generate offsprings...");
        let offsprings: Vec<Vec<usize>> = simulations
            .par_iter_mut()
            .zip(host_maps.par_iter())
            .map(|(simulation, host_map)| {
                simulation.mutate_infectants(host_map);
                simulation.replicate_infectants(host_map)
            })
            .collect();

        // transfer between compartments
        log::debug!("Transfer between compartments...");
        let transfer = plan.get_transfer_matrix(generation);
        let populations: Vec<Population> = (0..args.n_compartments)
            .into_par_iter()
            .map(|target| {
                (0..args.n_compartments)
                    .map(|origin| {
                        simulations[origin].subsample_population(
                            &offsprings[origin],
                            *transfer.get(target, origin),
                        )
                    })
                    .collect()
            })
            .collect();

        // update populations
        log::debug!("Update populations...");
        for (simulation, population) in simulations.iter_mut().zip(populations) {
            simulation.set_population(population);
        }
    }

    if let Some(bar) = bar {
        bar.finish_with_message("Done.");
    }
    log::info!("Finished simulation.");
}

#[cfg(not(feature = "parallel"))]
fn run(args: &Args, simulations: &mut [BasicSimulation], plan: Plan) {
    let bar = match args.disable_progress_bar {
        true => None,
        false => {
            let bar = ProgressBar::new(args.generations as u64);
            bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "[{bar:40}] {pos:>7}/{len:7} [{elapsed_precise} / {duration_precise}] {msg}",
                )
                .expect("Unable to create template.")
                .progress_chars("=> "),
        );
            Some(bar)
        }
    };

    for generation in 0..=args.generations {
        // logging
        log::debug!("Generate logging message for generation {generation}...");
        let population_sizes: Vec<usize> = simulations
            .iter()
            .map(|sim| sim.get_population().len())
            .collect();

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );

        if let Some(bar) = bar.as_ref() {
            bar.set_position(generation.try_into().unwrap());
            bar.set_message(format!("{population_sizes:?}"));
        }

        // write to output when sampling
        log::debug!("Process sampling...");
        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            sample(simulations, sample_size, generation, args);
        }

        // abort on last generation after sampling
        if generation == args.generations {
            break;
        }

        // adjust settings if needed
        if let Some(settings) = plan.get_settings(generation) {
            log::info!("Adjusting settings to:\n{}", settings);
            simulations.iter_mut().for_each(|simulation| {
                simulation.set_settings(settings.clone());
            });
        }

        // simulate compartmentalized population
        log::debug!("Generate offsprings...");
        let offsprings: Vec<Vec<usize>> = simulations
            .iter_mut()
            .map(|simulation| {
                simulation.increment_generation();
                if simulation.get_population().is_empty() {
                    return Vec::new();
                }
                let host_map = simulation.get_host_map();
                simulation.mutate_infectants(&host_map);
                simulation.replicate_infectants(&host_map)
            })
            .collect();

        // transfer between compartments
        log::debug!("Transfer between compartments...");
        let transfer = plan.get_transfer_matrix(generation);

        let populations: Vec<Population> = (0..args.n_compartments)
            .map(|target| {
                Population::from_iter((0..args.n_compartments).map(|origin| {
                    simulations[origin]
                        .subsample_population(&offsprings[origin], *transfer.get(target, origin))
                }))
            })
            .collect();

        // update populations
        for (simulation, population) in simulations.iter_mut().zip(populations) {
            simulation.set_population(population);
        }
    }

    if let Some(bar) = bar {
        bar.finish_with_message("Done.");
    }
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

    let plan = read_transfer_plan(args.transfer_plan.as_str());
    let sequence = load_sequence(args.sequence.as_str());

    // create wildtype
    let wildtype = Wildtype::new(sequence.clone());

    // load simulation settings
    let settings = SimulationSettings::read_from_file(args.settings.as_str()).unwrap_or_else(
        |err| match err {
            SimulationSettingsError::IoError(err) => {
                eprintln!("Unable to open settings from file '{}'.", args.settings);
                eprintln!("Reason: {err}");
                std::process::exit(1);
            }
            SimulationSettingsError::YamlError(err) => {
                eprintln!("Unable to load settings from file '{}'.", args.settings);
                eprintln!("Reason: {err}");
                std::process::exit(1);
            }
        },
    );
    log::info!("Loaded settings\n{}", settings);

    // create and write fitness table
    let fitness_table = FitnessTable::from_model(&sequence, 4, settings.fitness_model.clone())
        .unwrap_or_else(|err| {
            eprintln!("Unable to create fitness table.");
            eprintln!("Reason: {}", err.message);
            std::process::exit(1);
        });
    let fitness_table2 = FitnessTable::from_model(&sequence, 4, settings.fitness_model.clone())
        .unwrap_or_else(|err| {
            eprintln!("Unable to create fitness table.");
            eprintln!("Reason: {}", err.message);
            std::process::exit(1);
        });
    let mut fitness_file =
        io::BufWriter::new(fs::File::create(args.fitness_table.clone()).unwrap());
    fitness_table.write(&mut fitness_file).unwrap();

    // create individual compartments
    println!("Creating {} compartments...", args.n_compartments);
    let mut simulations: Vec<BasicSimulation> =
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
