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
    population: Population,
    fitness_table: FitnessTable,
    simulation_settings: SimulationSettings,
    mutation_sampler: Binomial,
}

impl Simulation {
    pub fn new(
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
        for (_host, infectants) in host_map {
            for infectant in infectants {
                let n_mutations = self.mutation_sampler.sample(&mut rng) as usize;

                if n_mutations <= 0 {
                    continue;
                }

                let mut infectant_ref = Rc::clone(&self.population[*infectant]);
                let sequence_length = infectant_ref.borrow().get_length();
                let sites = (0..sequence_length).choose_multiple(&mut rng, n_mutations);
                for site in sites {
                    let base = infectant_ref.borrow().get_base(site);
                    match base {
                        Some(val) => {
                            let dist = WeightedIndex::new(
                                self.simulation_settings.substitution_matrix[val as usize],
                            )
                            .unwrap();
                            let new_base = dist.sample(&mut rng);
                            infectant_ref = Rc::clone(&infectant_ref)
                                .borrow_mut()
                                .create_descendant(site, new_base as u8);
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
                    Ok(dist) => dist.sample(&mut rand::thread_rng()),
                    // if fitness is 0 => no offspring
                    Err(_) => 0.,
                };
                offspring[*infectant] = (offspring_sample / n_infectants as f64) as usize;
            }
        }
        offspring
    }

    pub fn subsample_population(&mut self, offspring_map: &Vec<usize>) -> Population {
        let mut rng = rand::thread_rng();
        let offspring_size: usize = offspring_map.iter().sum();
        let sample_size = min(
            (offspring_size as f64 * self.simulation_settings.dilution) as usize,
            self.simulation_settings.max_population,
        );
        let sampler = match WeightedIndex::new(offspring_map) {
            Ok(s) => s,
            Err(_) => {
                println!("Population went extinct.");
                process::exit(1);
            }
        };
        (0..sample_size)
            .map(|_| Rc::clone(&self.population[sampler.sample(&mut rng)]))
            .collect()
    }

    pub fn next_generation(&mut self) {
        let infectant_map = self.get_infectant_map();
        let host_map = self.get_host_map(&infectant_map);
        self.mutate_infectants(&host_map);
        let offspring = self.replicate_infectants(&host_map);
        let population = self.subsample_population(&offspring);
        self.set_population(population);
    }

    pub fn print_population(&self) -> Vec<String> {
        self.population
            .iter()
            .map(|hap| hap.borrow().get_string())
            .collect::<Vec<String>>()
    }
}
