#![feature(let_chains)]
#![feature(once_cell)]
#![feature(test)]

extern crate test;

use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
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
use virolution::plan::*;
use virolution::population;
use virolution::population::Population;
use virolution::references::HaplotypeRef;
use virolution::simulation::*;
use virolution::simulation_settings::*;

fn setup(args: &Args) {
    // setup logger
    simple_logging::log_to_file(args.log_file.as_str(), log::LevelFilter::Info)
        .expect("Failed to init logging.");

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
    settings: SimulationSettings,
) -> Vec<Simulation> {
    (0..args.n_compartments)
        .into_iter()
        .map(|compartment_idx| {
            let init_population: Population = if compartment_idx == 0 {
                population![wildtype.clone(); args.initial_population_size]
            } else {
                Population::new()
            };
            Simulation::new(
                wildtype.get_clone(),
                init_population,
                fitness_table.clone(),
                settings.clone(),
                0,
            )
        })
        .collect()
}

fn sample(simulations: &[Simulation], sample_size: usize, generation: usize, args: &Args) {
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
            .into_iter()
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
        // logging
        let population_sizes: Vec<usize> = simulations
            .iter()
            .map(|sim| sim.get_population().len())
            .collect();

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );

        // update progress bar
        bar.set_position(generation.try_into().unwrap());
        bar.set_message(format!("{population_sizes:?}"));

        // write to output when sampling
        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            sample(simulations, sample_size, generation, args);
        }

        // abort on last generation after sampling
        if generation == args.generations {
            break;
        }

        // simulate compartmentalized population in parallel
        simulations.iter_mut().for_each(|simulation| {
            simulation.increment_generation();
        });
        let host_maps: Vec<HostMap> = simulations
            .par_iter()
            .map(|simulation| simulation.get_host_map())
            .collect();
        let offsprings: Vec<Vec<usize>> = simulations
            .par_iter_mut()
            .zip(host_maps.par_iter())
            .map(|(simulation, host_map)| {
                simulation.mutate_infectants(host_map);
                simulation.replicate_infectants(host_map)
            })
            .collect();

        // transfer between compartments
        let transfer = plan.get_transfer_matrix(generation);
        let populations: Vec<Population> = (0..args.n_compartments)
            .into_par_iter()
            .map(|target| {
                (0..args.n_compartments)
                    .map(|origin| {
                        simulations[origin]
                            .subsample_population(&offsprings[origin], transfer[target][origin])
                    })
                    .collect()
            })
            .collect();

        // update populations
        for (simulation, population) in simulations.iter_mut().zip(populations) {
            simulation.set_population(population);
        }
    }
    bar.finish_with_message("Done.");
    log::info!("Finished simulation.");
}

#[cfg(not(feature = "parallel"))]
fn run(args: &Args, simulations: &mut [Simulation], plan: Plan) {
    // init progress bar
    let bar = ProgressBar::new(args.generations as u64);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("[{bar:40}] {pos:>7}/{len:7} [{elapsed_precise} / {duration_precise}] {msg}")
            .expect("Unable to create template.")
            .progress_chars("=> "),
    );

    for generation in 0..=args.generations {
        // write to output when sampling
        let sample_size = plan.get_sample_size(generation);
        if sample_size > 0 {
            sample(simulations, sample_size, generation, args);
        }

        // adjust settings if needed
        plan.get_settings(generation).and_then(|settings| {
            for simulation in simulations.iter_mut() {
                simulation.set_settings(settings.clone());
            }
            Some(())
        });

        // simulate compartmentalized population
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
        let transfer = plan.get_transfer_matrix(generation);

        let populations: Vec<Population> = (0..args.n_compartments)
            .into_iter()
            .map(|target| {
                Population::from_iter((0..args.n_compartments).map(|origin| {
                    simulations[origin]
                        .subsample_population(&offsprings[origin], transfer[target][origin])
                }))
            })
            .collect();

        // update populations
        for (simulation, population) in simulations.iter_mut().zip(populations) {
            simulation.set_population(population);
        }

        // logging
        let population_sizes: Vec<usize> = simulations
            .iter()
            .map(|sim| sim.get_population().len())
            .collect();

        log::info!(
            r###"
    generation={generation}
    population_sizes={population_sizes:?}"###
        );

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
    use virolution::population;

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
        let fitness_model = FitnessModel::new(distribution, UtilityFunction::Linear);
        let fitness_table = FitnessTable::from_model(&sequence, 4, fitness_model.clone()).unwrap();

        let wt = Wildtype::new(sequence);
        let population: Population = population![wt.clone(); 10];
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
            fitness_model,
        };
        let mut simulation = Simulation::new(wt, population, fitness_table, settings, 0);
        b.iter(|| {
            simulation.next_generation();
        })
    }
}
