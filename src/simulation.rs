use super::fitness::FitnessTable;
use super::haplotype::{Haplotype, HaplotypeRef};
use rand::prelude::*;
use rand_distr::Bernoulli;
use std::cell::RefCell;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::rc::Rc;

type Population = Vec<HaplotypeRef>;

pub struct Simulation {
    population: Population,
    fitness_table: FitnessTable,
    host_population_size: usize,
    infection_fraction: f64,
}

impl Simulation {
    pub fn new(
        population: Population,
        fitness_table: FitnessTable,
        host_population_size: usize,
        infection_fraction: f64,
    ) -> Self {
        Self {
            population: population,
            fitness_table: fitness_table,
            host_population_size: host_population_size,
            infection_fraction: infection_fraction,
        }
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

    pub fn get_host_map(&self, infectant_map: &Vec<Option<usize>>) -> HashMap<usize, Vec<usize>> {
        let mut host_map: HashMap<usize, Vec<usize>> = HashMap::new();
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
}
