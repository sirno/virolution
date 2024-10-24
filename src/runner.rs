use anyhow::{anyhow, Result};

use indicatif::{ProgressBar, ProgressStyle};
use itertools::Itertools;
#[cfg(not(feature = "parallel"))]
use rand::prelude::*;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cell::RefCell;
use std::cmp::min;
use std::fs;
use std::ops::Range;
use std::path::Path;
use std::rc::Rc;

use crate::args::Args;
use crate::config::{FitnessModelField, Parameters, Settings};
use crate::core::haplotype::set_number_of_fitness_tables;
use crate::core::{Ancestry, FitnessProvider, Haplotype, Historian, Population};
use crate::encoding::Nucleotide as Nt;
use crate::encoding::Symbol;
use crate::readwrite::{CsvSampleWriter, FastaSampleWriter, SampleWriter};
use crate::readwrite::{HaplotypeIO, PopulationIO};
use crate::references::HaplotypeRef;
#[cfg(feature = "parallel")]
use crate::simulation::HostMap;
use crate::simulation::{BasicSimulation, SimulationTrait};
#[cfg(not(feature = "parallel"))]
use crate::stats::population::{PopulationDistance, PopulationFrequencies};

pub struct Runner {
    args: Args,
    settings: Settings,
    wildtype: HaplotypeRef<Nt>,
    simulations: Vec<Box<SimulationTrait<Nt>>>,
    ancestry: Option<Ancestry>,
}

impl Runner {
    pub fn new(args: Args) -> Result<Runner> {
        // create output directory
        if let Some(outdir) = args.outdir.as_ref() {
            fs::create_dir_all(outdir)?;
        }

        // setup logger and rayon when in parallel mode
        Self::setup_logger(&args);
        #[cfg(feature = "parallel")]
        Self::setup_rayon(&args);

        // load settings and reference
        let settings = Self::load_settings(&args.settings)?;
        let wildtype = Haplotype::load_wildtype(args.sequence.as_str())?;

        dbg!(&wildtype);

        // initialize fitness
        let settings_path = Path::new(&args.settings).parent();
        let provider_map =
            Self::create_fitness_providers(&settings, &wildtype.get_sequence(), settings_path)?;
        let providers = provider_map
            .iter()
            .map(|(_range, provider)| provider)
            .collect::<Vec<_>>();
        Self::write_fitness_tables(&providers, args.outdir.as_ref().map(Path::new));

        // perform sanity checks
        if !settings
            .schedule
            .check_transfer_table_sizes(args.n_compartments)
        {
            return Err(anyhow!(
                "Incompatible transfer tables (they might be too small)"
            ));
        }

        // create individual compartments
        println!("Creating {} compartments...", args.n_compartments);
        let simulations: Vec<Box<SimulationTrait<Nt>>> =
            Self::create_simulations(&args, &wildtype, &provider_map, &settings.parameters[0]);

        // check if ancestry is requested and supported
        if cfg!(feature = "parallel") && args.ancestry.is_some() {
            log::warn!("Ancestry is not supported in parallel mode.");
            std::process::exit(1);
        }
        // initialize ancestry if requested
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

        let sequence_path = self
            .args
            .outdir
            .as_ref()
            .map(Path::new)
            .unwrap_or(Path::new("./"));

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
        let log_path = args
            .outdir
            .as_ref()
            .map(|p| Path::new(p).join(args.logfile.as_str()))
            .unwrap_or_else(|| Path::new(args.logfile.as_str()).to_path_buf());

        simple_logging::log_to_file(log_path, log_level).unwrap_or_else(|_| {
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

    fn create_fitness_providers<S: Symbol>(
        settings: &Settings,
        sequence: &[S],
        path: Option<&Path>,
    ) -> Result<Vec<(Range<usize>, FitnessProvider<S>)>> {
        let default_settings = &settings.parameters[0];
        let fitness_tables = match &default_settings.fitness_model {
            FitnessModelField::SingleHost(fitness_model) => {
                let mut fitness_model = fitness_model.clone();

                // paths within fitness model should be relative to the settings file
                if let Some(path) = path {
                    fitness_model.prepend_path(path.to_str().unwrap());
                }

                let _ = set_number_of_fitness_tables(1).is_ok();
                vec![(
                    0..default_settings.host_population_size,
                    FitnessProvider::from_model(0, sequence, &fitness_model)?,
                )]
            }
            FitnessModelField::MultiHost(fitness_models) => {
                let mut fitness_tables = Vec::new();
                let mut lower = 0;

                let _ = set_number_of_fitness_tables(fitness_models.len()).is_ok();

                for (id, fitness_model_frac) in fitness_models.iter().enumerate() {
                    let mut fitness_model = fitness_model_frac.fitness_model.clone();

                    // paths within fitness model should be relative to the settings file
                    if let Some(path) = path {
                        fitness_model.prepend_path(path.to_str().unwrap());
                    }

                    let n_hosts = (fitness_model_frac.fraction
                        * default_settings.host_population_size as f64)
                        .round() as usize;
                    let upper = min(lower + n_hosts, settings.parameters[0].host_population_size);

                    fitness_tables.push((
                        lower..upper,
                        FitnessProvider::from_model(id, sequence, &fitness_model)?,
                    ));
                    lower = upper;
                }

                fitness_tables
            }
        };
        Ok(fitness_tables)
    }

    fn write_fitness_tables<S: Symbol>(providers: &[&FitnessProvider<S>], path: Option<&Path>) {
        let write_path = path.unwrap_or_else(|| Path::new("./"));

        // write fitness tables
        providers.iter().for_each(|provider| {
            provider.write(write_path).unwrap();
        });
    }

    fn create_simulations(
        args: &Args,
        wildtype: &HaplotypeRef<Nt>,
        fitness_providers: &[(Range<usize>, FitnessProvider<Nt>)],
        parameters: &Parameters,
    ) -> Vec<Box<SimulationTrait<Nt>>> {
        (0..args.n_compartments)
            .map(|_compartment_idx| {
                let init_population: Population<Nt> = {
                    let initial_population_size = match args.initial_population_size {
                        Some(size) => size,
                        None => parameters.max_population,
                    };
                    match &args.initial_population_file {
                        Some(file_name) => Population::read(file_name, wildtype.clone())
                            .unwrap_or_else(|e| {
                                eprintln!("Unable to read initial population file: {file_name}.");
                                eprintln!("Error: {e}.");
                                std::process::exit(1);
                            }),
                        None => population![wildtype.clone(); initial_population_size],
                    }
                };
                BasicSimulation::new(
                    wildtype.clone(),
                    init_population,
                    fitness_providers.to_vec(),
                    parameters.clone(),
                    0,
                )
            })
            .map(|sim| Box::new(sim) as Box<SimulationTrait<Nt>>)
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

        let samples_dir: String = self
            .args
            .outdir
            .as_ref()
            .map(|p| Path::new(p).join("samples").to_string_lossy().to_string())
            .unwrap_or_else(|| "samples".to_string());

        let sample_writer: Box<dyn SampleWriter<Nt>> = match self.args.sampling_format.as_str() {
            "fasta" => Box::new(
                FastaSampleWriter::new::<Nt>(
                    &self.args.name,
                    &samples_dir,
                    Some(historian.clone()),
                )
                .unwrap_or_else(|err| {
                    eprintln!("Unable to create sample writer: {err}.");
                    std::process::exit(1);
                }),
            ),
            "csv" => Box::new(
                CsvSampleWriter::new::<Nt>(&self.args.name, &samples_dir, Some(historian.clone()))
                    .unwrap_or_else(|err| {
                        eprintln!("Unable to create sample writer: {err}.");
                        std::process::exit(1);
                    }),
            ),
            &_ => {
                eprintln!("Unknown sampling format.");
                std::process::exit(1);
            }
        };

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
            let populations: Vec<Population<Nt>> = (0..self.args.n_compartments)
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

        let samples_dir: String = self
            .args
            .outdir
            .as_ref()
            .map(|p| Path::new(p).join("samples").to_string_lossy().to_string())
            .unwrap_or_else(|| "samples".to_string());

        let sample_writer: Box<dyn SampleWriter<Nt>> = match self.args.sampling_format.as_str() {
            "fasta" => Box::new(
                FastaSampleWriter::new::<Nt>(
                    &self.args.name,
                    &samples_dir,
                    Some(historian.clone()),
                )
                .unwrap_or_else(|err| {
                    eprintln!("Unable to create sample writer: {err}.");
                    std::process::exit(1);
                }),
            ),
            "csv" => Box::new(
                CsvSampleWriter::new::<Nt>(&self.args.name, &samples_dir, Some(historian.clone()))
                    .unwrap_or_else(|err| {
                        eprintln!("Unable to create sample writer: {err}.");
                        std::process::exit(1);
                    }),
            ),
            &_ => {
                eprintln!("Unknown sampling format.");
                std::process::exit(1);
            }
        };

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
            let populations: Vec<Population<Nt>> = (0..self.args.n_compartments)
                .map(|target| {
                    Population::from_iter((0..self.args.n_compartments).map(|origin| {
                        self.simulations[origin]
                            .get_population()
                            .select(&indices[target][origin])
                    }))
                })
                .collect();

            // perform stat collection
            if let Some(stats) = self.settings.schedule.get_event_value("stats", generation) {
                for stat in stats.split(',') {
                    match stat {
                        "diversity" => {
                            let frequencies: Vec<Vec<f64>> = self
                                .simulations
                                .iter()
                                .map(|sim| sim.get_population().frequencies())
                                .collect();
                            log::info!("frequencies={frequencies:?}");
                        }
                        "distance" => {
                            let distances: Vec<f64> = self
                                .simulations
                                .iter()
                                .cartesian_product(self.simulations.iter())
                                .map(|(sim1, sim2)| {
                                    sim1.get_population().distance(sim2.get_population())
                                })
                                .collect();
                            log::info!("distance={distances:?}");
                        }
                        _ => {}
                    }
                }
            }

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
