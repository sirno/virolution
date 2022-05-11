use super::fitness::FitnessTable;
use super::haplotype::HaplotypeRef;
use super::simulation_settings::SimulationSettings;
use rand::prelude::*;
use rand_distr::{Bernoulli, Binomial, Poisson, WeightedIndex};
use std::cmp::min;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::process;
use std::rc::Rc;

pub type Population = Vec<HaplotypeRef>;
pub type HostMap = HashMap<usize, Vec<usize>>;

pub struct Simulation {
    wildtype: HaplotypeRef,
    population: Population,
    fitness_table: FitnessTable,
    simulation_settings: SimulationSettings,
    mutation_sampler: Binomial,
}

impl Simulation {
    pub fn new(
        wildtype: HaplotypeRef,
        population: Population,
        fitness_table: FitnessTable,
        simulation_settings: SimulationSettings,
    ) -> Self {
        let mutation_sampler = Binomial::new(
            population[0].borrow().get_length() as u64,
            simulation_settings.mutation_rate,
        )
        .unwrap();
        Self {
            wildtype: wildtype,
            population: population,
            fitness_table: fitness_table,
            simulation_settings: simulation_settings,
            mutation_sampler: mutation_sampler,
        }
    }

    pub fn get_population(&self) -> Population {
        self.population.clone()
    }

    pub fn set_population(&mut self, population: Population) {
        self.population = population;
    }

    pub fn get_infectant_map(&self) -> Vec<Option<usize>> {
        let mut rng = rand::thread_rng();
        let infection_distribution =
            Bernoulli::new(self.simulation_settings.infection_fraction).unwrap();
        let infectant_map: Vec<Option<usize>> = (0..self.population.len())
            .map(|_| {
                if infection_distribution.sample(&mut rng) {
                    Some(rng.gen_range(0..self.simulation_settings.host_population_size))
                } else {
                    None
                }
            })
            .collect();
        infectant_map
    }

    pub fn get_host_map(&self, infectant_map: &Vec<Option<usize>>) -> HostMap {
        let mut host_map: HostMap = HashMap::new();
        for (infectant, host) in infectant_map.iter().enumerate() {
            if let Some(host_id) = *host {
                match host_map.entry(host_id) {
                    Entry::Vacant(e) => drop(e.insert(vec![infectant])),
                    Entry::Occupied(mut e) => drop(e.get_mut().push(infectant)),
                }
            }
        }
        host_map
    }

    pub fn mutate_infectants(&mut self, host_map: &HostMap) {
        // mutate infectants based on host cell assignment
        let mut rng = rand::thread_rng();
        let sequence_length = self.wildtype.borrow().get_length();
        let site_vector: Vec<usize> = (0..sequence_length).collect();
        let site_options: &[usize] = site_vector.as_slice();
        for (_host, infectants) in host_map {
            for infectant in infectants {
                let n_mutations = self.mutation_sampler.sample(&mut rng) as usize;

                if n_mutations <= 0 {
                    continue;
                }

                let mut infectant_ref = self.population[*infectant].get_clone();
                let sites = site_options.choose_multiple(&mut rng, n_mutations);
                for site in sites {
                    let base = infectant_ref.borrow().get_base(*site);
                    match base {
                        Some(val) => {
                            let dist = WeightedIndex::new(
                                self.simulation_settings.substitution_matrix[val as usize],
                            )
                            .unwrap();
                            let new_base = dist.sample(&mut rng);
                            infectant_ref = infectant_ref
                                .get_clone()
                                .borrow_mut()
                                .create_descendant(*site, new_base as u8);

                            if infectant_ref.borrow().get_fitness(&self.fitness_table) <= 0. {
                                continue;
                            }
                        }
                        None => {}
                    }
                }

                self.population[*infectant] = infectant_ref;
            }
        }
    }

    pub fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize> {
        // replicate infectants within each host
        let mut rng = rand::thread_rng();
        let mut offspring = vec![0; self.population.len()];
        for (_host, infectants) in host_map {
            let n_infectants = infectants.len();
            for infectant in infectants {
                let fitness = self.population[*infectant]
                    .borrow()
                    .get_fitness(&self.fitness_table);
                let offspring_sample = match Poisson::new(
                    fitness * self.simulation_settings.basic_reproductive_number,
                ) {
                    Ok(dist) => dist.sample(&mut rng),
                    // if fitness is 0 => no offspring
                    Err(_) => 0.,
                };
                offspring[*infectant] = (offspring_sample / n_infectants as f64) as usize;
            }
        }
        offspring
    }

    pub fn subsample_population(&mut self, offspring_map: &Vec<usize>, factor: f64) -> Population {
        let mut rng = rand::thread_rng();
        let offspring_size: usize = offspring_map.iter().sum();
        let sample_size = (factor
            * min(
                (offspring_size as f64 * self.simulation_settings.dilution) as usize,
                self.simulation_settings.max_population,
            ) as f64) as usize;
        let sampler = match WeightedIndex::new(offspring_map) {
            Ok(s) => s,
            Err(_) => {
                println!("Population went extinct.");
                process::exit(1);
            }
        };
        (0..sample_size)
            .map(|_| self.population[sampler.sample(&mut rng)].get_clone())
            .collect()
    }

    pub fn next_generation(&mut self) {
        let infectant_map = self.get_infectant_map();
        let host_map = self.get_host_map(&infectant_map);
        self.mutate_infectants(&host_map);
        let offspring = self.replicate_infectants(&host_map);
        let population = self.subsample_population(&offspring, 1.);
        self.set_population(population);
    }

    pub fn print_population(&self) -> Vec<String> {
        self.population
            .iter()
            .map(|hap| hap.borrow().get_string())
            .collect::<Vec<String>>()
    }
}
