//! Population module
//!
//! The population module contains the `Population` struct, which is a collection
//! of haplotypes. The `Population` struct is used to represent the population
//! while only storing a single copy of each haplotype. This is done by storing
//! a vector of indices that point to the haplotypes in the population.
//!

use itertools::Itertools;
use rand::prelude::*;
use rand_distr::WeightedAliasIndex;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::collections::HashMap;
use std::ops::Index;

use crate::references::HaplotypeRef;

pub type Haplotypes = HashMap<usize, HaplotypeRef>;

#[macro_export]
macro_rules! population {
    () => {
        $crate::core::Population::new()
    };
    ($haplotype:expr, $size:expr) => {
        $crate::core::Population::from_haplotype($haplotype, $size)
    };
    ($( $haplotype:expr ),+) => {
        {
        let mut population = $crate::core::Population::new();
        $(
            population.push($haplotype);
        )+
        population
    }
    };
    ($( $haplotype:expr ; $size:expr ),+) => {
        $crate::core::Population::from_iter(vec![$( $crate::core::Population::from_haplotype($haplotype, $size) ),+])
    };
}

/// A `Population` is a collection of haplotypes.
#[derive(Clone, Debug, Default, PartialEq)]
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
            .unwrap_or_else(|| panic!("No haplotype with index {index}"))
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
            .or_else(|| panic!("No haplotype with index {ref_id}"))
    }
}

impl FromIterator<Population> for Population {
    fn from_iter<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        let mut population: Vec<usize> = Vec::new();
        let mut haplotypes: Haplotypes = Haplotypes::new();

        for pop in iter {
            population.extend(pop.population);
            haplotypes.extend(pop.haplotypes);
        }

        Self {
            population,
            haplotypes,
        }
    }
}

impl Population {
    /// Construct a new, empty `Population`.
    pub fn new() -> Self {
        Self {
            population: Vec::new(),
            haplotypes: Haplotypes::new(),
        }
    }

    /// Construct a new `Population` from a `HaplotypeRef` and a size.
    pub fn from_haplotype(haplotype: HaplotypeRef, size: usize) -> Self {
        let ref_id = haplotype.get_id();
        let population = vec![ref_id; size];
        let mut haplotypes = Haplotypes::new();
        haplotypes.insert(haplotype.get_id(), haplotype);
        Self {
            population,
            haplotypes,
        }
    }

    /// Construct a new `Population` from a `Vec` of `HaplotypeRef`s.
    pub fn from_references(references: Vec<HaplotypeRef>) -> Self {
        let mut population = Vec::new();
        let mut haplotypes = Haplotypes::new();

        for reference in references {
            let ref_id = reference.get_id();
            population.push(ref_id);
            haplotypes.insert(ref_id, reference);
        }

        Self {
            population,
            haplotypes,
        }
    }

    /// Get an iterator over the `Population`.
    pub fn iter(&self) -> PopulationIterator {
        PopulationIterator {
            population: self,
            index: 0,
        }
    }

    /// Check if the `Population` is empty.
    pub fn is_empty(&self) -> bool {
        self.population.is_empty()
    }

    /// Get the size of the `Population`.
    pub fn len(&self) -> usize {
        self.population.len()
    }

    /// Insert a `HaplotypeRef` at a specific position.
    ///
    /// This will overwrite the haplotype at that position, but not remove its
    /// reference from the `Haplotypes` map.
    pub fn insert(&mut self, position: &usize, haplotype: &HaplotypeRef) {
        let ref_id = haplotype.get_id();
        self.population[*position] = ref_id;
        self.haplotypes
            .entry(ref_id)
            .or_insert_with(|| haplotype.clone());
    }

    /// Push a `HaplotypeRef` to the end of the `Population`.
    pub fn push(&mut self, haplotype: &HaplotypeRef) {
        let ref_id = haplotype.get_id();
        self.population.push(ref_id);
        self.haplotypes
            .entry(ref_id)
            .or_insert_with(|| haplotype.clone());
    }

    /// Get a reference to the `Population` vector.
    pub fn get_population(&self) -> &Vec<usize> {
        &self.population
    }

    /// Get a reference to the `Haplotypes` map.
    pub fn get_haplotypes(&self) -> &Haplotypes {
        &self.haplotypes
    }

    /// Choose `Population` of size `size` from the `Population`.
    pub fn choose_multiple(&self, rng: &mut ThreadRng, size: usize) -> Self {
        let population: Vec<usize> = self
            .population
            .choose_multiple(rng, size)
            .map(|&id| id)
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

    #[allow(dead_code)]
    /// Sanitize the `Population`.
    fn sanitize(&mut self) {
        self.haplotypes = self
            .population
            .iter()
            .unique()
            .map(|&id| (id, self.haplotypes[&id].clone()))
            .collect();
    }

    #[cfg(feature = "parallel")]
    /// Get a random subsample `Population` of the `Population`.
    pub fn sample(&self, size: usize, weights: &[usize]) -> Self {
        let sampler = WeightedAliasIndex::new(weights.to_vec()).unwrap();

        let population: Vec<usize> = (0..size)
            .into_par_iter()
            .map_init(rand::thread_rng, |rng, _| {
                self.population[sampler.sample(rng)]
            })
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

    #[cfg(not(feature = "parallel"))]
    /// Get a random subsample `Population` of the `Population`.
    pub fn sample(&self, size: usize, weights: &[usize]) -> Self {
        let sampler = WeightedAliasIndex::new(weights.to_vec()).unwrap();

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

    use crate::core::haplotype::Wildtype;

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
        population.insert(&0, &wt2);
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

    #[test]
    fn from_iter() {
        let wt1 = Wildtype::new(vec![Some(0x00); 10]);
        let mut population1 = Population::new();
        population1.push(&wt1);

        let wt2 = Wildtype::new(vec![Some(0x00); 10]);
        let mut population2 = Population::new();
        population2.push(&wt2);

        let population = Population::from_iter(vec![population1, population2]);
        assert_eq!(population.len(), 2);
        assert_eq!(population[&0], wt1);
        assert_eq!(population[&1], wt2);
    }

    #[test]
    fn macro_empty() {
        let population = population![];
        assert_eq!(population.len(), 0);
        assert!(population.is_empty());
    }

    #[test]
    fn macro_from_haplotype() {
        let wt1 = Wildtype::new(vec![Some(0x00); 10]);
        assert_eq!(population![wt1.clone(); 1].len(), 1);
        assert_eq!(population![wt1.clone(); 2].len(), 2);
        assert_eq!(population![wt1.clone(); 10].len(), 10);
        assert_eq!(population![wt1.clone(); 1_000_000].len(), 1_000_000);
    }

    #[test]
    fn macro_from_haplotypes() {
        let wt1 = Wildtype::new(vec![Some(0x00); 10]);
        let wt2 = Wildtype::new(vec![Some(0x00); 10]);
        let wt3 = Wildtype::new(vec![Some(0x00); 10]);
        let population = population![&wt1, &wt2, &wt3];

        assert_eq!(population.len(), 3);
        assert_eq!(population[&0], wt1);
        assert_eq!(population[&1], wt2);
        assert_eq!(population[&2], wt3);
    }

    #[test]
    fn macro_from_haplotypes_variable_lengths() {
        let wt1 = Wildtype::new(vec![Some(0x00); 10]);
        let wt2 = Wildtype::new(vec![Some(0x00); 10]);
        let wt3 = Wildtype::new(vec![Some(0x00); 10]);
        let population = population![wt1.clone(); 2, wt2.clone(); 1, wt3.clone(); 3];

        assert_eq!(population.len(), 6);
        assert_eq!(population[&0], wt1);
        assert_eq!(population[&1], wt1);
        assert_eq!(population[&2], wt2);
        assert_eq!(population[&3], wt3);
        assert_eq!(population[&4], wt3);
        assert_eq!(population[&5], wt3);
    }
}
