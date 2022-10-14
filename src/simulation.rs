use crate::fitness::FitnessTable;
use crate::haplotype::Haplotype;
use crate::references::sync::HaplotypeRef;
use crate::simulation_settings::SimulationSettings;
use itertools::Itertools;
use rand::prelude::*;
use rand_distr::{Bernoulli, Binomial, Poisson, WeightedIndex};
use rayon::prelude::*;
use std::cmp::min;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::sync::mpsc::channel;

pub type Population = Vec<HaplotypeRef>;
pub type HostMap = HashMap<usize, Vec<usize>>;

pub struct Simulation {
    wildtype: HaplotypeRef,
    population: Population,
    fitness_table: FitnessTable,
    simulation_settings: SimulationSettings,
    mutation_sampler: Binomial,
    recombination_sampler: Bernoulli,
    infection_sampler: Bernoulli,
}

impl Simulation {
    pub fn new(
        wildtype: HaplotypeRef,
        population: Population,
        fitness_table: FitnessTable,
        simulation_settings: SimulationSettings,
    ) -> Self {
        let mutation_sampler = Binomial::new(
            wildtype.borrow().get_length() as u64,
            simulation_settings.mutation_rate,
        )
        .unwrap();
        let recombination_sampler = Bernoulli::new(simulation_settings.recombination_rate).unwrap();
        let infection_sampler = Bernoulli::new(simulation_settings.infection_fraction).unwrap();
        Self {
            wildtype,
            population,
            fitness_table,
            simulation_settings,
            mutation_sampler,
            recombination_sampler,
            infection_sampler,
        }
    }

    pub fn get_population(&self) -> Population {
        self.population.clone()
    }

    pub fn set_population(&mut self, population: Population) {
        self.population = population;
    }

    pub fn get_infectant_map(&self) -> Vec<Option<usize>> {
        let infectant_map: Vec<Option<usize>> = (0..self.population.len())
            .into_par_iter()
            .map_init(
                || rand::thread_rng(),
                |rng, _| {
                    if self.infection_sampler.sample(rng) {
                        Some(rng.gen_range(0..self.simulation_settings.host_population_size))
                    } else {
                        None
                    }
                },
            )
            .collect();
        infectant_map
    }

    pub fn get_host_map(&self, infectant_map: &[Option<usize>]) -> HostMap {
        let mut host_map: HostMap = HashMap::new();
        for (infectant, host) in infectant_map.iter().enumerate() {
            if let Some(host_id) = *host {
                match host_map.entry(host_id) {
                    Entry::Vacant(e) => {
                        e.insert(vec![infectant]);
                    }
                    Entry::Occupied(mut e) => {
                        e.get_mut().push(infectant);
                    }
                }
            }
        }
        host_map
    }

    pub fn mutate_infectants(&mut self, host_map: &HostMap) {
        // mutate infectants based on host cell assignment
        let sequence_length = self.wildtype.borrow().get_length();
        let site_vector: Vec<usize> = (0..sequence_length).collect();
        let site_options: &[usize] = site_vector.as_slice();

        // create recombination channel
        let (recombination_sender, recombination_receiver) = channel();

        // recombine infectants
        host_map
            .into_par_iter()
            .for_each_with(recombination_sender, |sender, entry| {
                let mut rng = rand::thread_rng();

                let infectants = entry.1;
                let n_infectants = infectants.len();

                // Recombine
                if n_infectants > 1 {
                    for infectant_pair in infectants.iter().combinations(2) {
                        if self.recombination_sampler.sample(&mut rng) {
                            let mut pair = [infectant_pair[0], infectant_pair[1]];
                            pair.shuffle(&mut rng);
                            let infectant_a = pair[0];
                            let infectant_b = pair[1];
                            let mut positions =
                                rand::seq::index::sample(&mut rng, sequence_length, 2).into_vec();
                            positions.sort();
                            let recombinant = Haplotype::create_recombinant(
                                &self.population[*infectant_a],
                                &self.population[*infectant_b],
                                positions[0],
                                positions[1],
                            );
                            sender
                                .send((*positions.choose(&mut rng).unwrap(), recombinant))
                                .unwrap();
                        }
                    }
                }
            });

        // collect recominants
        for (position, recombinant) in recombination_receiver.iter() {
            self.population[position] = recombinant;
        }

        // create mutation channel
        let (mutation_sender, mutation_receiver) = channel();

        // mutate infectants
        host_map
            .into_par_iter()
            .for_each_with(mutation_sender, |sender, entry| {
                let mut rng = rand::thread_rng();

                let infectants = entry.1;

                for infectant in infectants {
                    let n_mutations = self.mutation_sampler.sample(&mut rng) as usize;

                    if n_mutations == 0 {
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

                    sender.send((*infectant, infectant_ref)).unwrap();
                }
            });

        // collect mutants
        for (position, mutant) in mutation_receiver.iter() {
            self.population[position] = mutant;
        }
    }

    pub fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize> {
        // replicate infectants within each host
        let (replicate_sender, replicate_receiver) = channel();

        host_map
            .into_par_iter()
            .for_each_with(replicate_sender, |sender, entry| {
                let mut rng = rand::thread_rng();

                let infectants = entry.1;
                let n_infectants = infectants.len() as f64;

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

                    sender
                        .send((*infectant, offspring_sample / n_infectants))
                        .unwrap();
                }
            });

        let mut offspring = vec![0; self.population.len()];
        for (infectant, offspring_sample) in replicate_receiver.iter() {
            offspring[infectant] = offspring_sample as usize;
        }
        offspring
    }

    pub fn subsample_population(&self, offspring_map: &Vec<usize>, factor: f64) -> Population {
        // if there is no offspring, return empty population
        if offspring_map.is_empty() {
            return Vec::new();
        }

        let mut rng = rand::thread_rng();
        let offspring_size: usize = offspring_map.iter().sum();
        let sample_size = (factor
            * min(
                (offspring_size as f64 * self.simulation_settings.dilution) as usize,
                self.simulation_settings.max_population,
            ) as f64) as usize;
        let sampler = WeightedIndex::new(offspring_map).unwrap();

        (0..sample_size)
            .map(|_| self.population[sampler.sample(&mut rng)].get_clone())
            .collect()
    }

    pub fn next_generation(&mut self) {
        // do nothing if there is no population
        if self.population.is_empty() {
            return;
        }

        // simulate infection
        let infectant_map = self.get_infectant_map();
        let host_map = self.get_host_map(&infectant_map);

        // simulate replication and mutation
        self.mutate_infectants(&host_map);
        let offspring = self.replicate_infectants(&host_map);

        // subsample population
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fitness::{ExponentialParameters, FitnessDistribution, MutationCategoryWeights};
    use crate::haplotype::Wildtype;

    const SIMULATION_SETTINGS: SimulationSettings = SimulationSettings {
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
    };

    const DISTRIBUTION: FitnessDistribution =
        FitnessDistribution::Exponential(ExponentialParameters {
            weights: MutationCategoryWeights {
                beneficial: 0.29,
                deleterious: 0.51,
                lethal: 0.2,
                neutral: 0.,
            },
            lambda_beneficial: 0.03,
            lambda_deleterious: 0.21,
        });

    #[test]
    fn next_generation() {
        let sequence = vec![Some(0x00); 100];

        let fitness_table = FitnessTable::new(&sequence, &4, DISTRIBUTION);

        let wt = Wildtype::new(sequence);
        let init_population: Population = (0..10).map(|_| wt.get_clone()).collect();
        let mut simulation =
            Simulation::new(wt, init_population, fitness_table, SIMULATION_SETTINGS);
        simulation.next_generation()
    }

    #[test]
    fn next_generation_without_population() {
        let sequence = vec![Some(0x00); 100];

        let fitness_table = FitnessTable::new(&sequence, &4, DISTRIBUTION);

        let wt = Wildtype::new(sequence);
        let init_population: Population = (0..0).map(|_| wt.get_clone()).collect();
        let mut simulation =
            Simulation::new(wt, init_population, fitness_table, SIMULATION_SETTINGS);

        simulation.next_generation()
    }
}
