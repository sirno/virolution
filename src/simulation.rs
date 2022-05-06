use super::fitness::FitnessTable;
use super::haplotype::{Haplotype, HaplotypeRef};
use rand::prelude::*;
use rand_distr::{Bernoulli, Binomial, Poisson, WeightedIndex};
use std::cmp::min;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::rc::Rc;

type Population = Vec<HaplotypeRef>;
type HostMap = HashMap<usize, Vec<usize>>;

pub struct Simulation {
    population: Population,
    fitness_table: FitnessTable,
    substitution_matrix: [[f64; 4]; 4],
    host_population_size: usize,
    infection_fraction: f64,
    basic_reproductive_number: f64,
    mutation_sampler: Binomial,
    max_population: usize,
    dilution: f64,
}

impl Simulation {
    pub fn new(
        population: Population,
        fitness_table: FitnessTable,
        host_population_size: usize,
        infection_fraction: f64,
        basic_reproductive_number: f64,
        mutation_rate: f64,
    ) -> Self {
        let mutation_sampler =
            Binomial::new(population[0].borrow().get_length() as u64, mutation_rate).unwrap();
        Self {
            population: population,
            fitness_table: fitness_table,
            substitution_matrix: [
                [0., 1., 1., 1.],
                [1., 0., 1., 1.],
                [1., 1., 0., 1.],
                [1., 1., 1., 0.],
            ],
            host_population_size: host_population_size,
            infection_fraction: infection_fraction,
            basic_reproductive_number: basic_reproductive_number,
            mutation_sampler: mutation_sampler,
            max_population: 100,
            dilution: 0.01,
        }
    }

    pub fn get_population(&self) -> Population {
        self.population.clone()
    }

    pub fn get_infectant_map(&self) -> Vec<Option<usize>> {
        let mut rng = rand::thread_rng();
        let infection_distribution = Bernoulli::new(self.infection_fraction).unwrap();
        let infectant_map: Vec<Option<usize>> = (0..self.population.len())
            .map(|_| {
                if infection_distribution.sample(&mut rng) {
                    Some(rng.gen_range(0..self.host_population_size))
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
                            let dist =
                                WeightedIndex::new(self.substitution_matrix[val as usize]).unwrap();
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
                let offspring_sample = match Poisson::new(fitness * self.basic_reproductive_number)
                {
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
        println!("offspring_size: {:?}", offspring_size);
        let sample_size = min(
            (offspring_size as f64 * self.dilution) as usize,
            self.max_population,
        );
        let sampler = WeightedIndex::new(offspring_map).unwrap();
        (0..sample_size)
            .map(|_| Rc::clone(&self.population[sampler.sample(&mut rng)]))
            .collect()
    }
}
