use super::haplotype::Haplotype;
use rand::prelude::*;
use std::cell::RefCell;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::rc::Rc;

type Population = Vec<Rc<RefCell<Haplotype>>>;

pub struct Simulation {
    population: Population,
    host_population_size: usize,
}

impl Simulation {
    fn get_infectant_map(&self) -> Vec<usize> {
        let mut rng = rand::thread_rng();
        let infectant_map: Vec<usize> = (0..self.population.len())
            .map(|_| rng.gen_range(0..self.host_population_size))
            .collect();
        infectant_map
    }

    fn get_host_map(&self, infectant_map: &Vec<usize>) -> HashMap<usize, Vec<usize>> {
        let mut host_map: HashMap<usize, Vec<usize>> = HashMap::new();
        for (infectant, host) in infectant_map.iter().enumerate() {
            match host_map.entry(*host) {
                Entry::Vacant(e) => drop(e.insert(vec![infectant])),
                Entry::Occupied(mut e) => drop(e.get_mut().push(infectant)),
            }
        }
        host_map
    }
}
