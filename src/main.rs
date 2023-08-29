#![feature(let_chains)]
#![feature(get_mut_unchecked)]
#![feature(extract_if)]
#![feature(test)]

extern crate test;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use seq_io::fasta;
use seq_io::fasta::Record;
use std::cell::RefCell;
use std::cmp::min;
use std::fs;
use std::io;
use std::ops::Range;
use std::panic::catch_unwind;
use std::rc::Rc;
use virolution::args::*;
use virolution::config::{FitnessModelField, Parameters, Schedule, Settings, SettingsError};
use virolution::fitness::*;
use virolution::haplotype::*;
use virolution::historian::Historian;
use virolution::population;
use virolution::population::Population;
use virolution::references::HaplotypeRef;
use virolution::sample_writer::*;
use virolution::simulation::*;

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
    fitness_tables: &[(Range<usize>, FitnessTable)],
    parameters: &Parameters,
) -> Vec<Box<SimulationTrait>> {
    (0..args.n_compartments)
        .map(|compartment_idx| {
            let init_population: Population = if compartment_idx == 0 {
                let initial_population_size = match args.initial_population_size {
                    Some(size) => size,
                    None => parameters.max_population,
                };
                population![wildtype.clone(); initial_population_size]
            } else {
                Population::new()
            };
            BasicSimulation::new(
                wildtype.clone(),
                init_population,
                fitness_tables.to_vec(),
                parameters.clone(),
                0,
            )
        })
        .map(|sim| Box::new(sim) as Box<SimulationTrait>)
        .collect()
}

#[cfg(feature = "parallel")]
fn run(args: &Args, simulations: &mut Vec<Box<SimulationTrait>>, schedule: Schedule) {
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

    let historian = Rc::new(RefCell::new(Historian::new()));
    unsafe {
        // Leak this memory to avoid a free after use
        let ptr = Rc::into_raw(historian.clone());
        Rc::increment_strong_count(ptr);
        Rc::from_raw(ptr);
    }

    let sample_writer: Box<dyn SampleWriter> = Box::new(
        FastaSampleWriter::new(&args.name, &args.outdir, Some(historian)).unwrap_or_else(|err| {
            eprintln!("Unable to create sample writer: {err}.");
            std::process::exit(1);
        }),
    );

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
        let sample_size = schedule.get_sample_size(generation);
        if sample_size > 0 {
            sample_writer
                .sample_from_simulations(simulations, sample_size)
                .unwrap_or_else(|err| {
                    eprintln!("Unable to write sample: {err}.");
                    std::process::exit(1);
                })
        }

        // abort on last generation after sampling
        if generation == args.generations {
            break;
        }

        // adjust settings if needed
        if let Some(parameters) = schedule.get_settings(generation) {
            log::info!("Adjusting settings to:\n{}", settings);
            simulations.iter_mut().for_each(|simulation| {
                simulation.set_parameters(parameters.clone());
            });
        }

        // increment generation
        simulations.iter_mut().for_each(|simulation| {
            simulation.increment_generation();
        });

        // simulate compartmentalized population in parallel
        log::debug!("Generate host maps...");
        let host_maps: Vec<HostMap> = simulations
            .par_iter_mut()
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
        let transfer = schedule.get_transfer_matrix(generation);
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
fn run(args: &Args, simulations: &mut [Box<SimulationTrait>], schedule: Schedule) {
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

    let historian = Rc::new(RefCell::new(Historian::new()));
    unsafe {
        // Leak this memory to avoid a free after use
        let ptr = Rc::into_raw(historian.clone());
        Rc::increment_strong_count(ptr);
        Rc::from_raw(ptr);
    }

    let sample_writer: Box<dyn SampleWriter> = Box::new(
        FastaSampleWriter::new(&args.name, &args.outdir, Some(historian.clone())).unwrap_or_else(
            |err| {
                eprintln!("Unable to create sample writer: {err}.");
                std::process::exit(1);
            },
        ),
    );

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
        let sample_size = schedule.get_sample_size(generation);
        if sample_size > 0 {
            sample_writer
                .sample_from_simulations(simulations, sample_size)
                .unwrap_or_else(|err| {
                    eprintln!("Unable to write sample: {err}.");
                    std::process::exit(1);
                });
        }

        // abort on last generation after sampling
        if generation == args.generations {
            break;
        }

        // adjust settings if needed
        if let Some(settings) = schedule.get_settings(generation) {
            log::info!("Adjusting settings to:\n{}", settings);
            simulations.iter_mut().for_each(|simulation| {
                simulation.set_parameters(settings.clone());
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
        let transfer = schedule.get_transfer_matrix(generation);

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

    // load simulation settings
    let settings: Settings =
        Settings::read_from_file(args.settings.as_str()).unwrap_or_else(|err| match err {
            SettingsError::IoError(err) => {
                eprintln!("Unable to open settings from file '{}'.", args.settings);
                eprintln!("Reason: {err}");
                std::process::exit(1);
            }
            SettingsError::YamlError(err) => {
                eprintln!("Unable to load settings from file '{}'.", args.settings);
                eprintln!("Reason: {err}");
                std::process::exit(1);
            }
        });
    log::info!("Loaded settings\n{}", settings);

    // load sequence
    let sequence = load_sequence(args.sequence.as_str());

    // create wildtype
    let wildtype = Wildtype::new(sequence.clone());

    // create and write fitness table
    let fitness_tables = match &settings.parameters[0].fitness_model {
        FitnessModelField::SingleHost(fitness_model) => {
            vec![(
                0..settings.parameters[0].host_population_size,
                FitnessTable::from_model(&sequence, 4, fitness_model.clone()).unwrap_or_else(
                    |err| {
                        eprintln!("Unable to create fitness table.");
                        eprintln!("Reason: {}", err.message);
                        std::process::exit(1);
                    },
                ),
            )]
        }
        FitnessModelField::MultiHost(fitness_models) => {
            let mut fitness_tables = Vec::new();
            let mut lower = 0;
            for fitness_model_frac in fitness_models {
                let fitness_model = fitness_model_frac.fitness_model.clone();
                let n_hosts = (fitness_model_frac.fraction
                    * settings.parameters[0].host_population_size as f64)
                    .round() as usize;
                let upper = min(lower + n_hosts, settings.parameters[0].host_population_size);
                fitness_tables.push((
                    lower..upper,
                    FitnessTable::from_model(&sequence, 4, fitness_model).unwrap_or_else(|err| {
                        eprintln!("Unable to create fitness table.");
                        eprintln!("Reason: {}", err.message);
                        std::process::exit(1);
                    }),
                ));
                lower = upper;
            }
            fitness_tables
        }
    };

    for (idx, (_, fitness_table)) in fitness_tables.clone().iter().enumerate() {
        let name = format!("fitness_table_{}.npy", idx);
        let mut fitness_file = io::BufWriter::new(fs::File::create(name).unwrap());
        fitness_table.write(&mut fitness_file).unwrap();
    }

    // create individual compartments
    println!("Creating {} compartments...", args.n_compartments);
    let mut simulations: Vec<Box<SimulationTrait>> = create_simulations(
        &args,
        &wildtype,
        fitness_tables.as_slice(),
        &settings.parameters[0],
    );

    // run simulation
    run(&args, &mut simulations, settings.schedule);

    // store tree if specified.
    log::info!("Storing tree...");
    if let Some(tree_file) = args.trees {
        fs::write(tree_file, wildtype.get_tree())
            .unwrap_or_else(|_| eprintln!("Unable to write tree file."));
    }
    log::info!("Finished storing tree.");

    log::info!("Storing sequences...");
    for (compartment_id, compartment) in simulations.iter().enumerate() {
        let mut sequence_file = csv::WriterBuilder::new()
            .from_path(format!("final.{}.csv", compartment_id))
            .expect("Unable to create final sequence file.");

        sequence_file
            .write_record(["haplotype", "count"])
            .expect("Unable to write header to final sequence file.");

        for (haplotype_ref, haplotype_count) in compartment.get_population().iter().counts() {
            sequence_file
                .write_record(&[
                    haplotype_ref.as_ref().get_string(),
                    haplotype_count.to_string(),
                ])
                .expect("Unable to write to samples file.")
        }
    }
}
