use crate::fitness::FitnessTable;
use crate::haplotype::Haplotype;
use crate::population::Population;
use crate::references::HaplotypeRef;
use crate::simulation_settings::SimulationSettings;
use itertools::Itertools;
use rand::prelude::*;
use rand_distr::{Bernoulli, Binomial, Poisson, WeightedIndex};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::HashMap;
#[cfg(feature = "parallel")]
use std::sync::mpsc::channel;

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
            wildtype.get_length() as u64,
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

    pub fn get_population(&self) -> &Population {
        &self.population
    }

    pub fn set_population(&mut self, population: Population) {
        self.population = population;
    }

    #[cfg(feature = "parallel")]
    pub fn get_infectant_map(&self) -> Vec<Option<usize>> {
        (0..self.population.len())
            .into_par_iter()
            .map_init(rand::thread_rng, |rng, _| {
                if self.infection_sampler.sample(rng) {
                    Some(rng.gen_range(0..self.simulation_settings.host_population_size))
                } else {
                    None
                }
            })
            .collect()
    }

    #[cfg(not(feature = "parallel"))]
    pub fn get_infectant_map(&self) -> Vec<Option<usize>> {
        let mut rng = rand::thread_rng();
        let infectant_map: Vec<Option<usize>> = (0..self.population.len())
            .into_iter()
            .map(|_| {
                if self.infection_sampler.sample(&mut rng) {
                    Some(rng.gen_range(0..self.simulation_settings.host_population_size))
                } else {
                    None
                }
            })
            .collect();
        infectant_map
    }

    pub fn get_host_map(&self, infectant_map: &[Option<usize>]) -> HostMap {
        let mut host_map: HostMap = HashMap::new();
        let capacity = infectant_map.len() / self.simulation_settings.host_population_size + 1;
        for (infectant, host) in infectant_map.iter().enumerate() {
            if let Some(host_id) = *host {
                host_map
                    .entry(host_id)
                    .or_insert_with(|| Vec::with_capacity(capacity))
                    .push(infectant);
            }
        }
        host_map
    }

    fn _recombine_infectants(
        &self,
        sequence_length: usize,
        infectants: &Vec<usize>,
    ) -> Vec<(usize, HaplotypeRef)> {
        let n_infectants = infectants.len();

        if n_infectants <= 1 {
            return Vec::new();
        }

        // Recombine
        infectants
            .iter()
            .combinations(2)
            .filter(|_| self.recombination_sampler.sample(&mut rand::thread_rng()))
            .map(|infectant_pair| {
                let mut rng = rand::thread_rng();
                let mut recombination_sites =
                    rand::seq::index::sample(&mut rng, sequence_length, 2).into_vec();
                recombination_sites.sort();
                let recombinant = Haplotype::create_recombinant(
                    &self.population[infectant_pair[0]],
                    &self.population[infectant_pair[1]],
                    recombination_sites[0],
                    recombination_sites[1],
                );
                (**infectant_pair.choose(&mut rng).unwrap(), recombinant)
            })
            .collect()
    }

    fn _mutate_infectants(
        &self,
        sequence_length: usize,
        infectants: &[usize],
    ) -> Vec<(usize, HaplotypeRef)> {
        infectants
            .iter()
            .filter_map(|infectant| {
                let mut rng = rand::thread_rng();
                let infectant_ref = &self.population[infectant];
                let n_mutations = self.mutation_sampler.sample(&mut rng) as usize;
                if n_mutations == 0 {
                    return None;
                };
                let mut mutation_sites =
                    rand::seq::index::sample(&mut rng, sequence_length, n_mutations).into_vec();
                mutation_sites.sort();
                let bases = mutation_sites
                    .iter()
                    .map(|position| {
                        let mut rng = rand::thread_rng();
                        match infectant_ref.get_base(position) {
                            Some(base) => {
                                let dist = WeightedIndex::new(
                                    self.simulation_settings.substitution_matrix[base as usize],
                                )
                                .unwrap();
                                Some(dist.sample(&mut rng) as u8)
                            }
                            None => None,
                        }
                    })
                    .collect();

                let descendant = Haplotype::create_descendant(infectant_ref, mutation_sites, bases);
                Some((*infectant, descendant))
            })
            .collect()
    }

    #[cfg(feature = "parallel")]
    pub fn mutate_infectants(&mut self, host_map: &HostMap) {
        // mutate infectants based on host cell assignment
        let sequence_length = self.wildtype.get_length();

        if self.simulation_settings.recombination_rate > 0. {
            // create recombination channel
            let (recombination_sender, recombination_receiver) = channel();

            // recombine infectants
            host_map
                .par_iter()
                .for_each_with(recombination_sender, |sender, host| {
                    let infectants = host.1;

                    // Recombine
                    self._recombine_infectants(sequence_length, infectants)
                        .into_iter()
                        .for_each(|recombination| {
                            sender.send(recombination).unwrap();
                        });
                });

            // collect recominants
            for (position, recombinant) in recombination_receiver.iter() {
                self.population.insert(position, &recombinant);
            }
        }

        // create mutation channel
        let (mutation_sender, mutation_receiver) = channel();

        // mutate infectants
        host_map
            .par_iter()
            .for_each_with(mutation_sender, |sender, entry| {
                let infectants = entry.1;
                self._mutate_infectants(sequence_length, infectants)
                    .into_iter()
                    .for_each(|mutation| {
                        sender.send(mutation).unwrap();
                    });
            });

        // collect mutants
        for (position, mutant) in mutation_receiver.iter() {
            self.population.insert(position, &mutant);
        }
    }

    #[cfg(not(feature = "parallel"))]
    pub fn mutate_infectants(&mut self, host_map: &HostMap) {
        // mutate infectants based on host cell assignment
        let sequence_length = self.wildtype.get_length();

        if self.simulation_settings.recombination_rate > 0. {
            // recombine infectants
            for host in host_map {
                let infectants = host.1;
                self._recombine_infectants(sequence_length, infectants)
                    .into_iter()
                    .for_each(|(position, recombinant)| {
                        self.population.insert(position, &recombinant);
                    });
            }
        }

        // mutate infectants
        for host in host_map {
            let infectants = host.1;
            self._mutate_infectants(sequence_length, infectants)
                .into_iter()
                .for_each(|(position, mutant)| {
                    self.population.insert(position, &mutant);
                });
        }
    }

    fn _replicate_infectants(&self, infectants: &[usize]) -> Vec<(usize, f64)> {
        infectants
            .iter()
            .filter_map(|infectant| {
                let mut rng = rand::thread_rng();
                let fitness = self.population[infectant].get_fitness(&self.fitness_table);
                match Poisson::new(fitness * self.simulation_settings.basic_reproductive_number) {
                    Ok(dist) => Some((*infectant, dist.sample(&mut rng))),
                    // if fitness is 0 => no offspring
                    Err(_) => None,
                }
            })
            .collect()
    }

    #[cfg(feature = "parallel")]
    pub fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize> {
        // replicate infectants within each host
        let (replicate_sender, replicate_receiver) = channel();

        host_map
            .par_iter()
            .for_each_with(replicate_sender, |sender, host| {
                let infectants = host.1;
                self._replicate_infectants(infectants)
                    .into_iter()
                    .for_each(|replication| {
                        sender.send(replication).unwrap();
                    });
            });

        let mut offspring = vec![0; self.population.len()];
        for (infectant, offspring_sample) in replicate_receiver.iter() {
            offspring[infectant] = offspring_sample as usize;
        }
        offspring
    }

    #[cfg(not(feature = "parallel"))]
    pub fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize> {
        let mut offspring = vec![0; self.population.len()];
        host_map.iter().for_each(|host| {
            let infectants = host.1;
            self._replicate_infectants(infectants)
                .into_iter()
                .for_each(|replication| {
                    offspring[replication.0] = replication.1 as usize;
                });
        });
        offspring
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
        let offspring_map = self.replicate_infectants(&host_map);

        // subsample population
        self.population = self.population.subsample_population(
            &offspring_map,
            1.,
            self.simulation_settings.dilution,
            self.simulation_settings.max_population,
        );
    }

    pub fn subsample_population(&self, offspring_map: &[usize], factor: f64) -> Population {
        self.population.subsample_population(
            offspring_map,
            factor,
            self.simulation_settings.dilution,
            self.simulation_settings.max_population,
        )
    }

    pub fn print_population(&self) -> Vec<String> {
        self.population
            .iter()
            .map(|haplotype| haplotype.get_string())
            .collect::<Vec<String>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fitness::{
        ExponentialParameters, FitnessDistribution, FitnessModel, MutationCategoryWeights,
        UtilityFunction,
    };
    use crate::haplotype::Wildtype;

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

    const FITNESS_MODEL: FitnessModel = FitnessModel {
        distribution: DISTRIBUTION,
        utility: UtilityFunction::Linear,
    };

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
        fitness_model: FITNESS_MODEL,
    };

    #[test]
    fn next_generation() {
        let sequence = vec![Some(0x00); 100];

        let fitness_table = FitnessTable::new(&sequence, 4, FITNESS_MODEL);

        let wt = Wildtype::new(sequence);
        let population = Population::with_size(10, wt.get_clone());
        let mut simulation = Simulation::new(wt, population, fitness_table, SIMULATION_SETTINGS);
        simulation.next_generation()
    }

    #[test]
    fn next_generation_without_population() {
        let sequence = vec![Some(0x00); 100];

        let fitness_table = FitnessTable::new(&sequence, 4, FITNESS_MODEL);

        let wt = Wildtype::new(sequence);
        let mut population = Population::new();
        (0..10).for_each(|_| population.push(&wt));
        let mut simulation = Simulation::new(wt, population, fitness_table, SIMULATION_SETTINGS);

        simulation.next_generation()
    }
}
