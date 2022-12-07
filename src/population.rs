use itertools::Itertools;
use rand::prelude::*;
use rand_distr::WeightedIndex;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
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
        if self.index >= self.population.len() {
            return None;
        }

        let ref_id = self.population.population[self.index];
        self.index += 1;
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
        let mut population: Vec<usize> = Vec::new();
        let mut haplotypes: Haplotypes = Haplotypes::new();

        for pop in populations {
            population.extend(pop.population);
            haplotypes.extend(pop.haplotypes);
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
    pub fn sample(&self, offspring_map: &[usize], factor: f64, dilution: f64, max: usize) -> Self {
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
    pub fn sample(&self, size: usize, weights: &[usize]) -> Self {
        let sampler = WeightedIndex::new(weights).unwrap();

        let mut rng = rand::thread_rng();
        let population: Vec<usize> = (0..size)
            .map(|_| self.population[sampler.sample(&mut rng)])
            .collect();
        let haplotypes: Haplotypes = population
            .iter()
            .unique()
            .map(|&id| (id, self.haplotypes[&id].clone()))
            .collect();

        Self {
            population,
            haplotypes,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::haplotype::Wildtype;

    #[test]
    fn is_empty() {
        let mut population = Population::new();
        assert!(population.is_empty());

        let wt = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt);
        assert!(!population.is_empty());
    }

    #[test]
    fn len() {
        let mut population = Population::new();
        assert_eq!(population.len(), 0);

        let wt = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt);
        assert_eq!(population.len(), 1);
    }

    #[test]
    fn insert() {
        let mut population = Population::new();
        let wt = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt);
        assert_eq!(population.len(), 1);
        assert_eq!(population[&0], wt);

        let wt2 = Wildtype::new(vec![Some(0x00); 10]);
        population.insert(0, &wt2);
        assert_eq!(population.len(), 1);
        assert_ne!(population[&0], wt);
        assert_eq!(population[&0], wt2);
    }

    #[test]
    fn push() {
        let mut population = Population::new();
        let wt = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt);
        assert_eq!(population.len(), 1);
        assert_eq!(population[&0], wt);

        let wt2 = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt2);
        assert_eq!(population.len(), 2);
        assert_eq!(population[&1], wt2);
    }

    #[test]
    fn iterate() {
        let mut population = Population::new();
        let wt = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt);
        let wt2 = Wildtype::new(vec![Some(0x00); 10]);
        population.push(&wt2);

        let mut iter = population.iter();
        assert_eq!(iter.next().unwrap().get_id(), wt.get_id());
        assert_eq!(iter.next().unwrap().get_id(), wt2.get_id());
        assert_eq!(iter.next(), None);
    }
}
