use itertools::Itertools;
use rand::prelude::*;
use rand_distr::WeightedIndex;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashMap;
use std::ops::Index;

use crate::references::HaplotypeRef;

pub type Genotypes = HashMap<usize, HaplotypeRef>;

#[derive(Clone)]
pub struct Population {
    population: Vec<usize>,
    genotypes: Genotypes,
}

impl Index<&usize> for Population {
    type Output = HaplotypeRef;

    fn index(&self, index: &usize) -> &Self::Output {
        self.genotypes
            .get(index)
            .unwrap_or_else(|| panic!("No genotype with index {}", index))
    }
}

impl Population {
    pub fn new() -> Self {
        Self {
            population: Vec::new(),
            genotypes: Genotypes::new(),
        }
    }

    pub fn with_size(size: usize, genotype: HaplotypeRef) -> Self {
        let ref_id = genotype.get_id();
        let population = vec![ref_id; size];
        let mut genotypes = Genotypes::new();
        genotypes.insert(genotype.get_id(), genotype);
        Self {
            population,
            genotypes,
        }
    }

    pub fn from(population: Vec<usize>, ancestors: &[Genotypes]) -> Self {
        let genotypes = HashMap::from_iter(population.iter().unique().map(|&id| {
            for ancestor in ancestors {
                match ancestor.get(&id) {
                    Some(haplotype) => return (id, haplotype.clone()),
                    None => continue,
                };
            }
            panic!("No haplotype with id {}", id);
        }));
        Self {
            population,
            genotypes,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.population.is_empty()
    }

    pub fn insert(&mut self, position: usize, haplotype: &HaplotypeRef) {
        let ref_id = haplotype.get_id();
        self.population[position] = ref_id;
        self.genotypes
            .entry(ref_id)
            .or_insert_with(|| haplotype.clone());
    }

    pub fn push(&mut self, haplotype: &HaplotypeRef) {
        let ref_id = haplotype.get_id();
        self.population.push(ref_id);
        self.genotypes
            .entry(ref_id)
            .or_insert_with(|| haplotype.clone());
    }

    pub fn get_genotypes(&self) -> Genotypes {
        self.genotypes.clone()
    }

    #[cfg(feature = "parallel")]
    pub fn subsample_population(
        &self,
        offspring_map: &[usize],
        factor: f64,
        dilution: f64,
        max: usize,
    ) -> Self {
        // if there is no offspring, return empty population
        if offspring_map.is_empty() {
            return Self::new();
        }

        let offspring_size: usize = offspring_map.par_iter().sum();
        let sample_size =
            (factor * min((offspring_size as f64 * dilution) as usize, max) as f64) as usize;
        let sampler = WeightedIndex::new(offspring_map).unwrap();

        let population = (0..sample_size)
            .into_par_iter()
            .map_init(rand::thread_rng, |rng, _| {
                self.population[sampler.sample(rng)]
            })
            .collect();

        Self::from(population, &[self.genotypes])
    }

    #[cfg(not(feature = "parallel"))]
    pub fn subsample_population(
        &self,
        offspring_map: &Vec<usize>,
        factor: f64,
        dilution: f64,
        max: usize,
    ) -> Self {
        // if there is no offspring, return empty population
        if offspring_map.is_empty() {
            return Self::new();
        }

        let offspring_size: usize = offspring_map.iter().sum();
        let sample_size =
            (factor * min((offspring_size as f64 * dilution) as usize, max) as f64) as usize;
        let sampler = WeightedIndex::new(offspring_map).unwrap();

        let mut rng = rand::thread_rng();
        let population = (0..sample_size)
            .map(|_| self.population[sampler.sample(&mut rng)])
            .collect();

        Self::from(population, &[self.genotypes])
    }
}
