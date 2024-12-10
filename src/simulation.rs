//! Simulation module
extern crate test;

use itertools::Itertools;
use rand::prelude::*;
use rand_distr::{Bernoulli, Binomial, Poisson, WeightedAliasIndex, WeightedIndex};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
#[cfg(feature = "parallel")]
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::{cmp::min_by, ops::Range};

use crate::config::Parameters;
use crate::core::hosts::Host;
use crate::core::{Haplotype, Population};
use crate::encoding::Symbol;
use crate::providers::Generation;
use crate::references::HaplotypeRef;

pub type HostMap = Vec<Vec<usize>>;
// TODO: Review host spec abstraction
pub type HostSpec = (Range<usize>, &'static str);

struct BasicHost {
    parameters: Parameters,

    infection_fitness_provider: &'static str,
    replicative_fitness_provider: &'static str,

    infection_sampler: Bernoulli,
    mutation_sampler: Binomial,
    recombination_sampler: Bernoulli,
    replication_sampler: Poisson<f64>,
}

impl BasicHost {
    pub fn new(
        n_sites: usize,
        parameters: &Parameters,
        replicative_fitness_provider: &'static str,
        infection_fitness_provider: &'static str,
    ) -> Self {
        Self {
            parameters: parameters.clone(),

            infection_fitness_provider,
            replicative_fitness_provider,

            infection_sampler: Bernoulli::new(parameters.infection_fraction).unwrap(),
            mutation_sampler: Binomial::new(n_sites as u64, parameters.mutation_rate).unwrap(),
            recombination_sampler: Bernoulli::new(parameters.recombination_rate).unwrap(),
            replication_sampler: Poisson::new(parameters.basic_reproductive_number).unwrap(),
        }
    }
}

impl<S: Symbol> Host<S> for BasicHost {
    fn infect(&self, _haplotype: &HaplotypeRef<S>) -> bool {
        true
    }

    fn mutate(&self, haplotype: &mut [HaplotypeRef<S>]) {}

    fn replicate(&self, haplotypes: &[HaplotypeRef<S>], offspring: &mut [usize]) {
        // let fitness_values: Vec<f64> = haplotypes
        //     .iter()
        //     .map(|infectant| {
        //         f64::try_from(
        //             infectant
        //                 .get_attribute_or_compute(self.replicative_fitness_provider)
        //                 .unwrap(),
        //         )
        //         .unwrap()
        //     })
        //     .collect();

        // let fitness_sum: f64 = fitness_values.iter().sum();
        // host.iter()
        //     .zip(fitness_values.iter())
        //     .for_each(|(infectant, fitness)| {
        //         if let Ok(dist) =
        //             Poisson::new(fitness * self.parameters.basic_reproductive_number / fitness_sum)
        //         {
        //             offspring[*infectant] = dist.sample(&mut rng) as usize;
        //         }
        //     });
    }
}

#[cfg(feature = "parallel")]
pub type SimulationTrait<S> = dyn Simulation<S> + Send + Sync;
#[cfg(not(feature = "parallel"))]
pub type SimulationTrait<S> = dyn Simulation<S>;

pub trait Simulation<S: Symbol> {
    fn increment_generation(&mut self);

    fn set_parameters(&mut self, parameters: Parameters);

    fn get_generation(&self) -> usize;

    fn get_population(&self) -> &Population<S>;
    fn set_population(&mut self, population: Population<S>);

    fn get_host_map(&self) -> HostMap;
    fn mutate_infectants(&mut self, host_map: &HostMap);
    fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize>;

    // methods that require `offspring_map` for processing
    fn target_size(&self, offspring_map: &[usize]) -> f64;
    fn sample_indices(&self, offspring_map: &[usize], factor: usize) -> Vec<usize>;
    fn subsample_population(&self, offspring_map: &[usize], dilution: f64) -> Population<S>;

    fn print_population(&self) -> Vec<String>;

    fn next_generation(&mut self) {
        // increment generation counter
        self.increment_generation();

        // do nothing if there is no population
        if self.get_population().is_empty() {
            return;
        }

        // simulate infection
        let host_map = self.get_host_map();

        // simulate replication and mutation
        self.mutate_infectants(&host_map);
        let offspring_map = self.replicate_infectants(&host_map);

        // subsample population
        self.set_population(self.subsample_population(&offspring_map, 1.));
    }
}

pub struct BasicSimulation<S: Symbol> {
    wildtype: HaplotypeRef<S>,
    population: Population<S>,
    hosts: Vec<HostSpec>,
    parameters: Parameters,
    mutation_sampler: Binomial,
    recombination_sampler: Bernoulli,
    infection_sampler: Bernoulli,
    generation: Arc<Generation>,
}

impl<S: Symbol> BasicSimulation<S> {
    pub fn new(
        wildtype: HaplotypeRef<S>,
        population: Population<S>,
        hosts: Vec<HostSpec>,
        parameters: Parameters,
        generation: Arc<Generation>,
    ) -> Self {
        let mutation_sampler =
            Binomial::new(wildtype.get_length() as u64, parameters.mutation_rate).unwrap();
        let recombination_sampler = Bernoulli::new(parameters.recombination_rate).unwrap();
        let infection_sampler = Bernoulli::new(parameters.infection_fraction).unwrap();
        Self {
            wildtype,
            population,
            hosts,
            parameters,
            mutation_sampler,
            recombination_sampler,
            infection_sampler,
            generation,
        }
    }

    fn _recombine_infectants(
        &self,
        sequence_length: usize,
        infectants: &[usize],
    ) -> Vec<(usize, HaplotypeRef<S>)> {
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
    ) -> Vec<(usize, HaplotypeRef<S>)> {
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
                        let dist = WeightedIndex::new(
                            self.parameters.substitution_matrix
                                [*infectant_ref.get_base(position).index()],
                        )
                        .unwrap();
                        S::decode(&(dist.sample(&mut rng) as u8))
                    })
                    .collect();

                let descendant = Haplotype::create_descendant(infectant_ref, mutation_sites, bases);
                Some((*infectant, descendant))
            })
            .collect()
    }
}

impl<S: Symbol> Simulation<S> for BasicSimulation<S> {
    fn increment_generation(&mut self) {
        self.generation.increment();
    }

    fn get_generation(&self) -> usize {
        self.generation.get()
    }

    fn get_population(&self) -> &Population<S> {
        &self.population
    }

    fn set_population(&mut self, population: Population<S>) {
        self.population = population;
    }

    /// Set the parameters for the simulation
    ///
    /// WARNING: Currently this method does only reinitialize part of the simulation. Do not use
    /// without revising its implementation.
    fn set_parameters(&mut self, parameters: Parameters) {
        self.parameters = parameters;
        self.mutation_sampler = Binomial::new(
            self.wildtype.get_length() as u64,
            self.parameters.mutation_rate,
        )
        .unwrap();
        self.recombination_sampler = Bernoulli::new(self.parameters.recombination_rate).unwrap();
        self.infection_sampler = Bernoulli::new(self.parameters.infection_fraction).unwrap();
    }

    fn get_host_map(&self) -> HostMap {
        let capacity = self.population.len() / self.parameters.host_population_size + 1;
        let mut host_map: HostMap =
            vec![Vec::with_capacity(capacity); self.parameters.host_population_size];
        let mut rng = rand::thread_rng();
        (0..self.population.len()).for_each(|infectant| {
            if self.infection_sampler.sample(&mut rng) {
                let host_id = rng.gen_range(0..self.parameters.host_population_size);
                host_map[host_id].push(infectant);
            }
        });
        host_map
    }

    #[cfg(feature = "parallel")]
    fn mutate_infectants(&mut self, host_map: &HostMap) {
        // mutate infectants based on host cell assignment
        let sequence_length = self.wildtype.get_length();

        if self.parameters.recombination_rate > 0. {
            // create recombination channel
            let (recombination_sender, recombination_receiver) = channel();

            // recombine infectants
            host_map
                .par_iter()
                .for_each_with(recombination_sender, |sender, infectants| {
                    // Recombine
                    self._recombine_infectants(sequence_length, infectants)
                        .into_iter()
                        .for_each(|recombination| {
                            sender.send(recombination).unwrap();
                        });
                });

            // collect recombinants
            for (position, recombinant) in recombination_receiver.iter() {
                self.population.insert(&position, &recombinant);
            }
        }

        // create mutation channel
        let (mutation_sender, mutation_receiver) = channel();

        // mutate infectant
        host_map
            .par_iter()
            .for_each_with(mutation_sender, |sender, infectants| {
                self._mutate_infectants(sequence_length, infectants)
                    .into_iter()
                    .for_each(|mutation| {
                        sender.send(mutation).unwrap();
                    });
            });

        // collect mutants
        for (position, mutant) in mutation_receiver.iter() {
            self.population.insert(&position, &mutant);
        }
    }

    #[cfg(not(feature = "parallel"))]
    fn mutate_infectants(&mut self, host_map: &HostMap) {
        // mutate infectants based on host cell assignment
        let sequence_length = self.wildtype.get_length();

        if self.parameters.recombination_rate > 0. {
            // recombine infectants
            host_map.iter().for_each(|infectants: &Vec<usize>| {
                self._recombine_infectants(sequence_length, infectants)
                    .iter()
                    .for_each(|(position, recombinant)| {
                        self.population.insert(position, recombinant);
                    });
            });
        }

        // mutate infectants
        host_map.iter().for_each(|infectants: &Vec<usize>| {
            self._mutate_infectants(sequence_length, infectants)
                .iter()
                .for_each(|(position, mutant)| {
                    self.population.insert(position, mutant);
                });
        });
    }

    #[cfg(feature = "parallel")]
    fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize> {
        let mut offspring = vec![0; self.population.len()];
        let offspring_ptr = offspring.as_mut_ptr() as usize;
        self.hosts.par_iter().for_each(|(range, provider_id)| {
            host_map[range.clone()].par_iter().for_each(|host| {
                let fitness_values: Vec<f64> = host
                    .iter()
                    .map(|infectant| {
                        f64::try_from(
                            self.population[infectant]
                                .get_attribute_or_compute(provider_id)
                                .unwrap(),
                        )
                        .unwrap()
                    })
                    .collect();
                let fitness_sum: f64 = fitness_values.iter().sum();
                host.par_iter()
                    .zip(fitness_values.iter())
                    .for_each(|(infectant, fitness)| {
                        if let Ok(dist) = Poisson::new(
                            fitness * self.parameters.basic_reproductive_number / fitness_sum,
                        ) {
                            unsafe {
                                (offspring_ptr as *mut usize)
                                    .add(*infectant)
                                    .write(dist.sample(&mut rand::thread_rng()) as usize);
                            }
                        }
                    });
            });
        });
        offspring
    }

    #[cfg(not(feature = "parallel"))]
    fn replicate_infectants(&self, host_map: &HostMap) -> Vec<usize> {
        let mut offspring = vec![0; self.population.len()];
        let mut rng = rand::thread_rng();
        self.hosts.iter().for_each(|(range, provider_id)| {
            host_map[range.clone()].iter().for_each(|host| {
                let fitness_values: Vec<f64> = host
                    .iter()
                    .map(|infectant| {
                        f64::try_from(
                            self.population[infectant]
                                .get_attribute_or_compute(provider_id)
                                .unwrap(),
                        )
                        .unwrap()
                    })
                    .collect();
                let fitness_sum: f64 = fitness_values.iter().sum();
                host.iter()
                    .zip(fitness_values.iter())
                    .for_each(|(infectant, fitness)| {
                        if let Ok(dist) = Poisson::new(
                            fitness * self.parameters.basic_reproductive_number / fitness_sum,
                        ) {
                            offspring[*infectant] = dist.sample(&mut rng) as usize;
                        }
                    });
            });
        });
        offspring
    }

    /// Compute the target size for offspring dilution
    ///
    /// The return value is f64 so callers can decide how to round the values before
    /// subsampling
    fn target_size(&self, offspring_map: &[usize]) -> f64 {
        if offspring_map.is_empty() {
            return 0.0;
        }

        let offspring_size: usize = offspring_map.iter().sum();

        min_by(
            offspring_size as f64 * self.parameters.dilution,
            self.parameters.max_population as f64,
            |a, b| a.partial_cmp(b).unwrap(),
        )
    }

    fn sample_indices(&self, offspring_map: &[usize], amount: usize) -> Vec<usize> {
        if offspring_map.is_empty() {
            return Vec::new();
        }

        let mut rng = rand::thread_rng();
        let sampler = WeightedAliasIndex::new(offspring_map.to_vec()).unwrap();

        (0..amount).map(|_| sampler.sample(&mut rng)).collect()
    }

    fn subsample_population(&self, offspring_map: &[usize], factor: f64) -> Population<S> {
        // if there is no offspring, return empty population
        if offspring_map.is_empty() {
            return Population::new();
        }

        let sample_size = (factor * self.target_size(offspring_map)) as usize;

        log::info!("Subsampling population to size {}", sample_size);

        self.population.sample(sample_size, offspring_map)
    }

    fn print_population(&self) -> Vec<String> {
        self.population
            .iter()
            .map(|haplotype| haplotype.get_string())
            .collect::<Vec<String>>()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::config::FitnessModelField;
    use crate::core::attributes::AttributeSetDefinition;
    use crate::core::fitness::init::{
        ExponentialParameters, FitnessDistribution, FitnessModel, MutationCategoryWeights,
    };
    use crate::core::fitness::utility::UtilityFunction;
    use crate::core::fitness::FitnessProvider;
    use crate::core::haplotype::Wildtype;
    use crate::encoding::Nucleotide as Nt;
    use std::sync::Arc;
    use test::Bencher;

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

    const POPULATION_SIZE: usize = 1000;

    const SETTINGS: Parameters = Parameters {
        mutation_rate: 1e-3,
        recombination_rate: 1e-5,
        substitution_matrix: [
            [0., 1., 1., 1.],
            [1., 0., 1., 1.],
            [1., 1., 0., 1.],
            [1., 1., 1., 0.],
        ],
        host_population_size: 5,
        infection_fraction: 1.0,
        basic_reproductive_number: 100.,
        max_population: POPULATION_SIZE,
        dilution: 0.17,
        fitness_model: FitnessModelField::SingleHost(FITNESS_MODEL),
    };

    fn setup_test_simulation() -> BasicSimulation<Nt> {
        let sequence = vec![Nt::A; 100];

        let mut attribute_definitions = AttributeSetDefinition::new();
        let name = "fitness";
        attribute_definitions.register(Arc::new(
            FitnessProvider::from_model(name, &sequence, &FITNESS_MODEL).unwrap(),
        ));
        let hosts = vec![(0..SETTINGS.host_population_size, name)];

        let wt = Wildtype::new(sequence, &attribute_definitions);
        let population: Population<Nt> = crate::population![wt.clone(), POPULATION_SIZE];
        let generation = Arc::new(Generation::new(0));
        BasicSimulation::new(wt, population, hosts, SETTINGS, generation)
    }

    #[test]
    fn increment_generation() {
        let mut simulation = setup_test_simulation();
        simulation.increment_generation();
        assert_eq!(simulation.generation.get(), 1);
    }

    #[test]
    fn get_population() {
        let simulation = setup_test_simulation();
        let population = simulation.get_population().clone();
        assert_eq!(simulation.population, population);
    }

    #[test]
    fn set_population() {
        let mut simulation = setup_test_simulation();
        let population: Population<Nt> = crate::population![simulation.wildtype.clone(), 42];
        simulation.set_population(population.clone());
        assert_eq!(simulation.population.len(), 42);
        assert_eq!(simulation.population, population);
    }

    #[test]
    fn get_host_map() {
        let simulation = setup_test_simulation();
        let host_map = simulation.get_host_map();
        assert_eq!(host_map.len(), SETTINGS.host_population_size);
        let n_infectants: usize = host_map.iter().map(|v| v.len()).sum();
        assert_eq!(n_infectants, SETTINGS.max_population);
    }

    #[test]
    fn next_generation() {
        let mut simulation = setup_test_simulation();
        simulation.next_generation()
    }

    #[test]
    fn next_generation_without_population() {
        let mut simulation = setup_test_simulation();
        simulation.set_population(Population::new());
        simulation.next_generation()
    }

    #[bench]
    fn bench_next_generation(b: &mut Bencher) {
        let mut simulation = setup_test_simulation();
        b.iter(|| simulation.next_generation());
    }

    #[bench]
    fn bench_get_host_map(b: &mut Bencher) {
        let simulation = setup_test_simulation();
        b.iter(|| simulation.get_host_map());
    }

    #[bench]
    fn bench_replicate_infectants(b: &mut Bencher) {
        let simulation = setup_test_simulation();
        let host_map = simulation.get_host_map();
        b.iter(|| simulation.replicate_infectants(&host_map));
    }

    #[bench]
    fn bench_subsample_population(b: &mut Bencher) {
        let simulation = setup_test_simulation();
        let host_map = simulation.get_host_map();
        let offspring_map = simulation.replicate_infectants(&host_map);
        b.iter(|| simulation.subsample_population(&offspring_map, 1.0));
    }
}
