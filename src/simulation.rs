//! Simulation module
extern crate test;

use itertools::Itertools;
use rand::prelude::*;
use rand_distr::weighted::{WeightedAliasIndex, WeightedIndex};
use rand_distr::{Bernoulli, Binomial, Poisson, Uniform};
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cmp::min_by;
use std::sync::Arc;

use crate::config::Parameters;
use crate::core::hosts::{Host, HostMapBuffer, HostSpecs};
use crate::core::population::Store;
use crate::core::{Haplotype, Population};
use crate::encoding::Symbol;
use crate::providers::Generation;
use crate::references::HaplotypeRef;

#[derive(Debug, Clone)]
pub struct BasicHost {
    parameters: Parameters,

    infection_fitness_provider: Option<&'static str>,
    replicative_fitness_provider: Option<&'static str>,

    site_sampler: Uniform<usize>,
    mutation_sampler: Binomial,
    recombination_sampler: Bernoulli,
    // replication_sampler: Poisson<f64>,
}

impl BasicHost {
    pub fn new(
        n_sites: usize,
        parameters: &Parameters,
        replicative_fitness_provider: Option<&'static str>,
        infection_fitness_provider: Option<&'static str>,
    ) -> Self {
        Self {
            parameters: parameters.clone(),

            infection_fitness_provider,
            replicative_fitness_provider,

            site_sampler: Uniform::new(0, n_sites).unwrap(),
            mutation_sampler: Binomial::new(n_sites as u64, parameters.mutation_rate).unwrap(),
            recombination_sampler: Bernoulli::new(parameters.recombination_rate).unwrap(),
            // replication_sampler: Poisson::new(parameters.basic_reproductive_number).unwrap(),
        }
    }
}

impl<S: Symbol> Host<S> for BasicHost {
    fn infect(&self, haplotype: &HaplotypeRef<S>, rng: &mut ThreadRng) -> bool {
        let infectivity = match self.infection_fitness_provider {
            Some(provider) => {
                let fitness = haplotype.get_or_compute_attribute(provider).unwrap();
                fitness.try_into().unwrap()
            }
            None => 1.,
        };
        rng.random::<f64>() < infectivity
    }

    fn mutate(&self, haplotype: &mut [HaplotypeRef<S>], rng: &mut ThreadRng) {
        if haplotype.len() > 1 && self.parameters.recombination_rate > 0. {
            // recombine infectants
            for (i, j) in (0..haplotype.len()).tuple_combinations() {
                if self.recombination_sampler.sample(rng) {
                    let mut recombination_sites =
                        [self.site_sampler.sample(rng), self.site_sampler.sample(rng)];
                    // ensure that recombination sites are different
                    while recombination_sites[0] == recombination_sites[1] {
                        recombination_sites[1] = self.site_sampler.sample(rng);
                    }
                    // sort recombination sites
                    recombination_sites.sort();
                    let recombinant = Haplotype::create_recombinant(
                        &haplotype[i],
                        &haplotype[j],
                        recombination_sites[0],
                        recombination_sites[1],
                    );
                    // randomly choose one of the parents to replace
                    haplotype[*[i, j].choose(rng).unwrap()] = recombinant;
                }
            }
        };

        // mutate infectants
        haplotype.iter_mut().for_each(|infectant| {
            let n_mutations = self.mutation_sampler.sample(rng) as usize;
            if n_mutations == 0 {
                return;
            };
            let mut mutation_sites = self
                .site_sampler
                .sample_iter(&mut *rng)
                .take(n_mutations)
                .collect_vec();
            mutation_sites.sort();
            let bases = mutation_sites
                .iter()
                .map(|position| {
                    let dist = WeightedIndex::new(
                        self.parameters.substitution_matrix[*infectant.get_base(position).index()],
                    )
                    .unwrap();
                    S::decode(&(dist.sample(rng) as u8))
                })
                .collect();

            let descendant = Haplotype::create_descendant(infectant, mutation_sites, bases);
            *infectant = descendant;
        });
    }

    fn replicate(
        &self,
        haplotypes: &[HaplotypeRef<S>],
        offspring: &mut [usize],
        rng: &mut ThreadRng,
    ) {
        assert_eq!(haplotypes.len(), offspring.len());

        // create float view to use as temporary fitness buffer
        let offspring_float_view = unsafe {
            std::slice::from_raw_parts_mut(offspring.as_mut_ptr() as *mut f64, offspring.len())
        };

        // collect fitness and compute sum
        for (infectant, offspring_field) in haplotypes.iter().zip(offspring_float_view.iter_mut()) {
            *offspring_field = match self.replicative_fitness_provider {
                Some(provider) => infectant
                    .get_or_compute_attribute(provider)
                    .unwrap()
                    .try_into()
                    .unwrap(),
                None => 1.,
            };
        }
        let fitness_sum: f64 = offspring_float_view.iter().sum();

        // sample offspring
        for i in 0..offspring.len() {
            if let Ok(dist) = Poisson::new(
                offspring_float_view[i] * self.parameters.basic_reproductive_number / fitness_sum,
            ) {
                offspring[i] = dist.sample(rng) as usize;
            }
        }
    }

    fn get_attributes(&self) -> Vec<&'static str> {
        let mut attributes = Vec::new();
        if let Some(attr) = self.infection_fitness_provider {
            attributes.push(attr);
        }
        if let Some(attr) = self.replicative_fitness_provider {
            attributes.push(attr);
        }
        attributes
    }

    fn clone_box(&self) -> Box<dyn Host<S>> {
        Box::new(self.clone())
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

    fn get_population(&self) -> &Population<Store<S>>;
    fn set_population(&mut self, population: Population<Store<S>>);

    fn get_host_specs(&self) -> &HostSpecs<S>;

    fn infect(&mut self);
    fn mutate_infectants(&mut self);
    fn replicate_infectants(&self) -> Vec<usize>;

    // methods that require `offspring_map` for processing
    fn target_size(&self, offspring_map: &[usize]) -> f64;
    fn sample_indices(&self, offspring_map: &[usize], factor: usize) -> Vec<usize>;
    fn subsample_population(&self, offspring_map: &[usize], dilution: f64) -> Population<Store<S>>;

    fn print_population(&self) -> Vec<String>;

    fn next_generation(&mut self) {
        // increment generation counter
        self.increment_generation();

        // do nothing if there is no population
        if self.get_population().is_empty() {
            return;
        }

        // simulate infection
        self.infect();

        // simulate replication and mutation
        self.mutate_infectants();
        let offspring_map = self.replicate_infectants();

        // subsample population
        self.set_population(self.subsample_population(&offspring_map, 1.));
    }
}

pub struct BasicSimulation<S: Symbol> {
    wildtype: HaplotypeRef<S>,
    population: Population<Store<S>>,
    parameters: Parameters,
    host_specs: HostSpecs<S>,
    mutation_sampler: Binomial,
    recombination_sampler: Bernoulli,
    generation: Arc<Generation>,

    // internal buffers
    host_map_buffer: HostMapBuffer,
}

impl<S: Symbol> BasicSimulation<S> {
    pub fn new(
        wildtype: HaplotypeRef<S>,
        population: Population<Store<S>>,
        parameters: Parameters,
        host_specs: HostSpecs<S>,
        generation: Arc<Generation>,
    ) -> Self {
        let mutation_sampler =
            Binomial::new(wildtype.get_length() as u64, parameters.mutation_rate).unwrap();
        let recombination_sampler = Bernoulli::new(parameters.recombination_rate).unwrap();
        let host_map_buffer = HostMapBuffer::new(parameters.host_population_size, population.len());
        Self {
            wildtype,
            population,
            parameters,
            host_specs,
            mutation_sampler,
            recombination_sampler,
            generation,
            host_map_buffer,
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
            .filter(|_| self.recombination_sampler.sample(&mut rand::rng()))
            .map(|infectant_pair| {
                let mut rng = rand::rng();
                let mut recombination_sites =
                    rand::seq::index::sample(&mut rng, sequence_length, 2).into_vec();
                recombination_sites.sort();
                let recombinant = Haplotype::create_recombinant(
                    &self.population.get(infectant_pair[0]),
                    &self.population.get(infectant_pair[1]),
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
                let mut rng = rand::rng();
                let infectant_ref = &self.population.get(infectant);
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
                        let mut rng = rand::rng();
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

    fn get_population(&self) -> &Population<Store<S>> {
        &self.population
    }

    fn set_population(&mut self, population: Population<Store<S>>) {
        self.population = population;

        // update host map buffer
        self.host_map_buffer.limit_infectants(self.population.len());
    }

    fn get_host_specs(&self) -> &HostSpecs<S> {
        &self.host_specs
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
    }

    fn infect(&mut self) {
        let host_sampler = Uniform::new(0, self.parameters.host_population_size)
            .expect("Invalid host population size");
        let mut rng = rand::rng();
        self.host_map_buffer
            .build(|(infectant, ref mut infection)| {
                let host_candidate = host_sampler.sample(&mut rng);
                **infection = self
                    .host_specs
                    .try_get_spec_from_index(host_candidate)
                    .and_then(|spec| {
                        if spec.host.infect(&self.population.get(&infectant), &mut rng) {
                            Some(host_candidate)
                        } else {
                            None
                        }
                    })
            });
    }

    #[cfg(feature = "parallel")]
    fn mutate_infectants(&mut self) {
        // mutate infectants based on host cell assignment
        self.host_specs.par_iter().for_each(|spec| {
            spec.range.clone().into_par_iter().for_each_init(
                || rand::rng(),
                |rng, position| {
                    let infectant_ids = self.host_map_buffer.get_slice(position);
                    let mut infectants = infectant_ids
                        .iter()
                        .map(|infectant_id| self.population.get(infectant_id).clone())
                        .collect::<Vec<HaplotypeRef<S>>>();
                    spec.host.mutate(infectants.as_mut_slice(), rng);
                    for (id, infectant) in infectant_ids.iter().zip(infectants.iter()) {
                        self.population.insert_sync(id, infectant.clone());
                    }
                },
            );
        });
    }

    #[cfg(not(feature = "parallel"))]
    fn mutate_infectants(&mut self) {
        // mutate infectants based on host cell assignment
        let mut rng = rand::rng();
        // recombine infectants
        self.host_specs.iter().for_each(|spec| {
            spec.range.clone().for_each(|position| {
                let infectant_ids = self.host_map_buffer.get_slice(position);
                let mut infectants = infectant_ids
                    .iter()
                    .map(|infectant_id| self.population.get(infectant_id).clone())
                    .collect::<Vec<HaplotypeRef<S>>>();
                spec.host.mutate(infectants.as_mut_slice(), &mut rng);
                for (id, infectant) in infectant_ids.iter().zip(infectants.iter()) {
                    self.population.insert(id, infectant.clone());
                }
            });
        });
    }

    #[cfg(feature = "parallel")]
    fn replicate_infectants(&self) -> Vec<usize> {
        // TODO: parallelize this method
        let mut offspring = vec![0; self.population.len()];
        let mut rng = rand::rng();
        self.host_specs.iter().for_each(|spec| {
            self.host_map_buffer
                .iter_range(spec.range.clone())
                .for_each(|infectants| {
                    let mut offspring_buf = vec![0; infectants.len()];
                    spec.host.replicate(
                        infectants
                            .iter()
                            .map(|infectant_id| self.population.get(infectant_id).clone())
                            .collect::<Vec<HaplotypeRef<S>>>()
                            .as_slice(),
                        offspring_buf.as_mut_slice(),
                        &mut rng,
                    );
                    // this can potentially be optimized by using the buffer directly and later
                    // reordering the values with the host_map.
                    for (infectant_id, n_offspring) in infectants.iter().zip(offspring_buf.iter()) {
                        offspring[*infectant_id] = *n_offspring;
                    }
                });
        });
        offspring
    }

    #[cfg(not(feature = "parallel"))]
    fn replicate_infectants(&self) -> Vec<usize> {
        let mut offspring = vec![0; self.population.len()];
        let mut rng = rand::rng();
        self.host_specs.iter().for_each(|spec| {
            self.host_map_buffer
                .iter_range(spec.range.clone())
                .for_each(|infectants| {
                    let mut offspring_buf = vec![0; infectants.len()];
                    spec.host.replicate(
                        infectants
                            .iter()
                            .map(|infectant_id| self.population.get(infectant_id).clone())
                            .collect::<Vec<HaplotypeRef<S>>>()
                            .as_slice(),
                        offspring_buf.as_mut_slice(),
                        &mut rng,
                    );
                    // this can potentially be optimized by using the buffer directly and later
                    // reordering the values with the host_map.
                    for (infectant_id, n_offspring) in infectants.iter().zip(offspring_buf.iter()) {
                        offspring[*infectant_id] = *n_offspring;
                    }
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

        let mut rng = rand::rng();
        let sampler = WeightedAliasIndex::new(offspring_map.to_vec()).unwrap();

        (0..amount).map(|_| sampler.sample(&mut rng)).collect()
    }

    fn subsample_population(&self, offspring_map: &[usize], factor: f64) -> Population<Store<S>> {
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
    use crate::config::{HostFitness, HostModel};
    use crate::core::fitness::utility::UtilityFunction;
    use crate::core::haplotype::Wildtype;
    use crate::encoding::Nucleotide as Nt;
    use crate::init::fitness::{
        ExponentialParameters, FitnessDistribution, FitnessModel, MutationCategoryWeights,
    };
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

    const FITNESS_MODEL: HostFitness = HostFitness {
        reproductive: Some(FitnessModel {
            distribution: DISTRIBUTION,
            utility: UtilityFunction::Linear,
        }),
        infective: None,
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
        host_population_size: 20,
        basic_reproductive_number: 100.,
        max_population: POPULATION_SIZE,
        dilution: 1.,
        host_model: HostModel::SingleHost(FITNESS_MODEL),
        n_hits: 1,
    };

    fn setup_test_simulation() -> BasicSimulation<Nt> {
        let sequence = vec![Nt::A; 100];

        let (attribute_definitions, host_specs) = SETTINGS
            .host_model
            .make_definitions(&SETTINGS, &sequence, None);

        let wt = Wildtype::new(sequence, &attribute_definitions);
        let population: Population<Store<Nt>> = crate::population![wt.clone(), POPULATION_SIZE];
        let generation = Arc::new(Generation::new(0));

        BasicSimulation::new(wt, population, SETTINGS, host_specs, generation)
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
        let population: Population<Store<Nt>> = crate::population![simulation.wildtype.clone(), 42];
        simulation.set_population(population.clone());
        assert_eq!(simulation.population.len(), 42);
        assert_eq!(simulation.population, population);
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
    fn bench_infect(b: &mut Bencher) {
        let mut simulation = setup_test_simulation();
        b.iter(|| simulation.infect());
    }

    #[bench]
    fn bench_replicate_infectants(b: &mut Bencher) {
        let mut simulation = setup_test_simulation();
        simulation.infect();
        b.iter(|| simulation.replicate_infectants());
    }

    #[bench]
    fn bench_subsample_population(b: &mut Bencher) {
        let mut simulation = setup_test_simulation();
        simulation.infect();
        let offspring_map = simulation.replicate_infectants();
        b.iter(|| simulation.subsample_population(&offspring_map, 1.0));
    }
}
