use itertools::Itertools;
use rand::prelude::*;
use rand_distr::WeightedIndex;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cmp::min;
use std::collections::HashMap;
use std::ops::Index;

use crate::references::HaplotypeRef;

pub type Haplotypes = HashMap<usize, HaplotypeRef>;

#[derive(Clone)]
pub struct Population {
    population: Vec<usize>,
    haplotypes: Haplotypes,
}

pub struct PopulationIterator<'a> {
    population: &'a Population,
    index: usize,
}

impl Index<&usize> for Population {
    type Output = HaplotypeRef;

    fn index(&self, index: &usize) -> &Self::Output {
        let ref_id = &self.population[*index];
        self.haplotypes
            .get(ref_id)
            .unwrap_or_else(|| panic!("No haplotype with index {}", index))
    }
}

impl<'a> Iterator for PopulationIterator<'a> {
    type Item = &'a HaplotypeRef;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index > self.population.len() {
            return None;
        }

        self.index += 1;
        let ref_id = self.population.population[self.index];
        self.population
            .haplotypes
            .get(&ref_id)
            .or_else(|| panic!("No haplotype with index {}", ref_id))
    }
}

impl Population {
    pub fn new() -> Self {
        Self {
            population: Vec::new(),
            haplotypes: Haplotypes::new(),
        }
    }

    pub fn with_size(size: usize, haplotype: HaplotypeRef) -> Self {
        let ref_id = haplotype.get_id();
        let population = vec![ref_id; size];
        let mut haplotypes = Haplotypes::new();
        haplotypes.insert(haplotype.get_id(), haplotype);
        Self {
            population,
            haplotypes,
        }
    }

    pub fn from_iter(populations: impl Iterator<Item = Self>) -> Self {
        let mut population = Vec::new();
        let mut haplotypes = Haplotypes::new();
        for pop in populations {
            pop.population.iter().unique().for_each(|&ref_id| {
                haplotypes.insert(ref_id, pop.haplotypes[&ref_id].clone());
            });
            population.extend(pop.population);
        }
        Self {
            population,
            haplotypes,
        }
    }

    pub fn iter(&self) -> PopulationIterator {
        PopulationIterator {
            population: self,
            index: 0,
        }
    }

    pub fn is_empty(&self) -> bool {
        self.population.is_empty()
    }

    pub fn len(&self) -> usize {
        self.population.len()
    }

    pub fn insert(&mut self, position: usize, haplotype: &HaplotypeRef) {
        let ref_id = haplotype.get_id();
        self.population[position] = ref_id;
        self.haplotypes
            .entry(ref_id)
            .or_insert_with(|| haplotype.clone());
    }

    pub fn push(&mut self, haplotype: &HaplotypeRef) {
        let ref_id = haplotype.get_id();
        self.population.push(ref_id);
        self.haplotypes
            .entry(ref_id)
            .or_insert_with(|| haplotype.clone());
    }

    pub fn get_haplotypes(&self) -> Haplotypes {
        self.haplotypes.clone()
    }

    pub fn choose_multiple(&self, rng: &mut ThreadRng, amount: usize) -> Vec<&HaplotypeRef> {
        self.population
            .choose_multiple(rng, amount)
            .map(|&id| {
                self.haplotypes
                    .get(&id)
                    .unwrap_or_else(|| panic!("No haplotype with index {}", id))
            })
            .collect()
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

        Self {
            population,
            haplotypes: self.haplotypes.clone(),
        }
    }

    #[cfg(not(feature = "parallel"))]
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

        let offspring_size: usize = offspring_map.iter().sum();
        let sample_size =
            (factor * min((offspring_size as f64 * dilution) as usize, max) as f64) as usize;
        let sampler = WeightedIndex::new(offspring_map).unwrap();

        let mut rng = rand::thread_rng();
        let population = (0..sample_size)
            .map(|_| self.population[sampler.sample(&mut rng)])
            .collect();

        Self {
            population,
            haplotypes: self.haplotypes.clone(),
        }
    }
}
