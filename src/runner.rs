use anyhow::Result;

use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
use rand::prelude::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cell::RefCell;
use std::cmp::min;
use std::fs;
use std::io;
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;

use crate::args::Args;
use crate::config::{FitnessModelField, Parameters, Settings};
use crate::core::haplotype::set_number_of_fitness_tables;
use crate::core::{Ancestry, FitnessTable, Haplotype, Historian, Population};
use crate::readwrite::{FastaSampleWriter, SampleWriter};
use crate::readwrite::{HaplotypeIO, PopulationIO};
use crate::references::HaplotypeRef;
#[cfg(feature = "parallel")]
use crate::simulation::HostMap;
use crate::simulation::{BasicSimulation, SimulationTrait};

pub struct Runner {
    args: Args,
    settings: Settings,
    wildtype: HaplotypeRef,
    simulations: Vec<Box<SimulationTrait>>,
    ancestry: Option<Ancestry>,
}

impl Runner {
    pub fn new(args: Args) -> Result<Runner> {
        Self::setup_logger(&args);
        #[cfg(feature = "parallel")]
        Self::setup_rayon(&args);

        let settings = Self::load_settings(&args.settings)?;
        let wildtype = Haplotype::load_wildtype(args.sequence.as_str())?;

        let fitness_tables = Self::create_fitness_tables(&settings, &wildtype.get_sequence())?;
        Self::write_fitness_tables(&fitness_tables, Path::new(args.outdir.as_str()).parent());

        // create individual compartments
        println!("Creating {} compartments...", args.n_compartments);
        let simulations: Vec<Box<SimulationTrait>> = Self::create_simulations(
            &args,
            &wildtype,
            fitness_tables.as_slice(),
            &settings.parameters[0],
        );

        let ancestry: Option<Ancestry> = args.ancestry.as_ref().map(|_| Ancestry::new());

        Ok(Self {
            args,
            settings,
            wildtype,
            simulations,
            ancestry,
        })
    }

    pub fn start(&mut self) {
        self.run();
        self.finish();
    }

    fn finish(&self) {
        // store tree if specified.
        log::info!("Storing tree...");
        if let Some(tree_file) = &self.args.trees {
            fs::write(tree_file, self.wildtype.get_tree())
                .unwrap_or_else(|_| eprintln!("Unable to write tree file."));
        }
        log::info!("Finished storing tree.");

        log::info!("Storing sequences...");

        let sequence_path = Path::new(self.args.outdir.as_str())
            .parent()
            .unwrap_or_else(|| Path::new("./"));
        for (compartment_id, compartment) in self.simulations.iter().enumerate() {
            log::info!("Storing sequences for compartment {}...", compartment_id);
            let mut sequence_file = csv::WriterBuilder::new()
                .from_path(sequence_path.join(format!("final.{}.csv", compartment_id)))
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

        if let Some(ancestry) = &self.ancestry {
            log::info!("Storing ancestry...");
            let ancestry_file = self.args.ancestry.as_ref().unwrap();
            let names: Vec<String> = self
                .simulations
                .iter()
                .enumerate()
                .flat_map(|(t, s)| {
                    s.get_population()
                        .iter()
                        .map(move |h| format!("{}_{}", h.get_block_id(), t))
                })
                .collect();
            fs::write(ancestry_file, ancestry.get_tree(names.as_slice()))
                .unwrap_or_else(|_| eprintln!("Unable to write ancestry file."));
        }
    }

    /// Setup logging level and file
    fn setup_logger(args: &Args) {
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
    }

    /// Setup rayon thread pool
    #[cfg(feature = "parallel")]
    fn setup_rayon(args: &Args) {
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

    /// Load settings from file
    fn load_settings(path: &str) -> Result<Settings> {
        let settings: Settings = Settings::read_from_file(path)?;
        log::info!("Loaded settings\n{}", settings);
        Ok(settings)
    }

    fn create_fitness_tables(
        settings: &Settings,
        sequence: &Vec<Option<u8>>,
    ) -> Result<Vec<(Range<usize>, FitnessTable)>> {
        let fitness_tables = match &settings.parameters[0].fitness_model {
            FitnessModelField::SingleHost(fitness_model) => {
                let _ = set_number_of_fitness_tables(1);
                vec![(
                    0..settings.parameters[0].host_population_size,
                    FitnessTable::from_model(0, sequence, 4, fitness_model.clone())?,
                )]
            }
            FitnessModelField::MultiHost(fitness_models) => {
                let mut fitness_tables = Vec::new();
                let mut lower = 0;
                let _ = set_number_of_fitness_tables(fitness_models.len());
                for (id, fitness_model_frac) in fitness_models.iter().enumerate() {
                    let fitness_model = fitness_model_frac.fitness_model.clone();
                    let n_hosts = (fitness_model_frac.fraction
                        * settings.parameters[0].host_population_size as f64)
                        .round() as usize;
                    let upper = min(lower + n_hosts, settings.parameters[0].host_population_size);
                    fitness_tables.push((
                        lower..upper,
                        FitnessTable::from_model(id, sequence, 4, fitness_model)?,
                    ));
                    lower = upper;
                }
                fitness_tables
            }
        };
        Ok(fitness_tables)
    }

    fn write_fitness_tables(fitness_tables: &[(Range<usize>, FitnessTable)], path: Option<&Path>) {
        let sequence_path = path.unwrap_or_else(|| Path::new("./"));
        for (idx, (_, fitness_table)) in fitness_tables.iter().enumerate() {
            let name = format!("fitness_table_{}.npy", idx);
            let mut fitness_file =
                io::BufWriter::new(fs::File::create(sequence_path.join(name)).unwrap());
            fitness_table.write(&mut fitness_file).unwrap();
        }
    }

    fn create_simulations(
        args: &Args,
        wildtype: &HaplotypeRef,
        fitness_tables: &[(Range<usize>, FitnessTable)],
        parameters: &Parameters,
    ) -> Vec<Box<SimulationTrait>> {
        (0..args.n_compartments)
            .map(|_compartment_idx| {
                let init_population: Population = {
                    let initial_population_size = match args.initial_population_size {
                        Some(size) => size,
                        None => parameters.max_population,
                    };
                    match &args.initial_population_file {
                        Some(file_name) => Population::read(file_name, wildtype.clone()).unwrap(),
                        None => population![wildtype.clone(); initial_population_size],
                    }
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
    fn run(&mut self) {
        let bar = match self.args.disable_progress_bar {
            true => None,
            false => {
                let bar = ProgressBar::new(self.args.generations as u64);
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
            FastaSampleWriter::new(&self.args.name, &self.args.outdir, Some(historian))
                .unwrap_or_else(|err| {
                    eprintln!("Unable to create sample writer: {err}.");
                    std::process::exit(1);
                }),
        );

        for generation in 0..=self.args.generations {
            // logging
            log::debug!("Generate logging message for generation {generation}...");
            let population_sizes: Vec<usize> = self
                .simulations
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
            let sample_size = self.settings.schedule.get_sample_size(generation);
            if sample_size > 0 {
                sample_writer
                    .sample_from_simulations(&self.simulations, sample_size)
                    .unwrap_or_else(|err| {
                        eprintln!("Unable to write sample: {err}.");
                        std::process::exit(1);
                    })
            }

            // abort on last generation after sampling
            if generation == self.args.generations {
                break;
            }

            // adjust settings if needed
            if let Some(parameters) = self.settings.schedule.get_settings(generation) {
                log::info!("Adjusting settings to:\n{}", parameters);
                self.simulations.iter_mut().for_each(|simulation| {
                    simulation.set_parameters(parameters.clone());
                });
            }

            // increment generation
            self.simulations.iter_mut().for_each(|simulation| {
                simulation.increment_generation();
            });

            // simulate compartmentalized population in parallel
            log::debug!("Generate host maps...");
            let host_maps: Vec<HostMap> = self
                .simulations
                .par_iter_mut()
                .map(|simulation| simulation.get_host_map())
                .collect();
            log::debug!("Generate offsprings...");
            let offsprings: Vec<Vec<usize>> = self
                .simulations
                .par_iter_mut()
                .zip(host_maps.par_iter())
                .map(|(simulation, host_map)| {
                    simulation.mutate_infectants(host_map);
                    simulation.replicate_infectants(host_map)
                })
                .collect();

            // transfer between compartments
            log::debug!("Transfer between compartments...");
            let transfer = self.settings.schedule.get_transfer_matrix(generation);
            let populations: Vec<Population> = (0..self.args.n_compartments)
                .into_par_iter()
                .map(|target| {
                    (0..self.args.n_compartments)
                        .map(|origin| {
                            self.simulations[origin].subsample_population(
                                &offsprings[origin],
                                *transfer.get(target, origin),
                            )
                        })
                        .collect()
                })
                .collect();

            // update populations
            log::debug!("Update populations...");
            for (simulation, population) in self.simulations.iter_mut().zip(populations) {
                simulation.set_population(population);
            }
        }

        if let Some(bar) = bar {
            bar.finish_with_message("Done.");
        }
        log::info!("Finished simulation.");
    }

    #[cfg(not(feature = "parallel"))]
    fn run(&mut self) {
        let bar = match self.args.disable_progress_bar {
            true => None,
            false => {
                let bar = ProgressBar::new(self.args.generations as u64);
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
            FastaSampleWriter::new(&self.args.name, &self.args.outdir, Some(historian.clone()))
                .unwrap_or_else(|err| {
                    eprintln!("Unable to create sample writer: {err}.");
                    std::process::exit(1);
                }),
        );

        for generation in 0..=self.args.generations {
            // logging
            log::debug!("Generate logging message for generation {generation}...");
            let population_sizes: Vec<usize> = self
                .simulations
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
            let sample_size = self.settings.schedule.get_sample_size(generation);
            if sample_size > 0 {
                sample_writer
                    .sample_from_simulations(&self.simulations, sample_size)
                    .unwrap_or_else(|err| {
                        eprintln!("Unable to write sample: {err}.");
                        std::process::exit(1);
                    });
            }

            // abort on last generation after sampling
            if generation == self.args.generations {
                break;
            }

            // adjust settings if needed
            if let Some(settings) = self.settings.schedule.get_settings(generation) {
                log::info!("Adjusting settings to:\n{}", settings);
                self.simulations.iter_mut().for_each(|simulation| {
                    simulation.set_parameters(settings.clone());
                });
            }

            // simulate compartmentalized population
            log::debug!("Generate offsprings...");
            let offsprings: Vec<Vec<usize>> = self
                .simulations
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
            let transfer = self.settings.schedule.get_transfer_matrix(generation);

            // choose transfer amounts
            let target_sizes: Vec<f64> = (0..self.args.n_compartments)
                .map(|target| self.simulations[target].target_size(&offsprings[target]))
                .collect();
            let migration_amount: Vec<Vec<usize>> = (0..self.args.n_compartments)
                .map(|target| {
                    (0..self.args.n_compartments)
                        .map(|origin| {
                            if target == origin {
                                return 0;
                            }

                            let mut rng = rand::thread_rng();

                            let migration_rate = transfer.get(target, origin);
                            let migration_amount = target_sizes[target] * migration_rate;
                            let residue =
                                if rng.gen::<f64>() < migration_amount - migration_amount.floor() {
                                    1
                                } else {
                                    0
                                };
                            migration_amount as usize + residue
                        })
                        .collect()
                })
                .collect();

            // sample indices to transfer
            let indices: Vec<Vec<Vec<usize>>> = (0..self.args.n_compartments)
                .map(|target| {
                    (0..self.args.n_compartments)
                        .map(|origin| {
                            let transfer_amount = if target == origin {
                                target_sizes[target] as usize
                                    - migration_amount[target].iter().sum::<usize>()
                            } else {
                                migration_amount[target][origin]
                            };
                            self.simulations[origin]
                                .sample_indices(&offsprings[origin], transfer_amount)
                        })
                        .collect()
                })
                .collect();

            // store ancestry from indices
            // TODO: ONLY WORKS WITH CONSTANT POPULATION SIZE
            if let Some(ancestry) = &mut self.ancestry {
                let ancestors: Vec<usize> = indices
                    .iter()
                    .flat_map(|t| {
                        t.iter()
                            .enumerate()
                            .flat_map(|(i, v)| {
                                let n = self.simulations[..i]
                                    .iter()
                                    .map(|s| s.get_population().len())
                                    .sum::<usize>();
                                v.iter().map(|a| n + (*a)).collect::<Vec<usize>>()
                            })
                            .collect::<Vec<_>>()
                    })
                    .collect();
                ancestry.add_ancestors(ancestors);
            }

            // perform the migration and create new populations
            let populations: Vec<Population> = (0..self.args.n_compartments)
                .map(|target| {
                    Population::from_iter((0..self.args.n_compartments).map(|origin| {
                        self.simulations[origin]
                            .get_population()
                            .select(&indices[target][origin])
                    }))
                })
                .collect();

            // update populations
            for (simulation, population) in self.simulations.iter_mut().zip(populations) {
                simulation.set_population(population);
            }
        }

        if let Some(bar) = bar {
            bar.finish_with_message("Done.");
        }
        log::info!("Finished simulation.");
    }
}
