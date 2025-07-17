//! Population module.
//!
//! The population module contains the `Population` struct, which is a collection
//! of haplotypes. The `Population` struct is used to represent the population
//! while only storing a single copy of each haplotype. This is done by storing
//! a vector of indices that point to the haplotypes in the population.

#[cfg(feature = "parallel")]
use dashmap::DashMap;
use itertools::Itertools;
use rand::prelude::*;
use rand_distr::weighted::WeightedAliasIndex;
#[cfg(feature = "parallel")]
use rayon::prelude::*;
use std::cell::Cell;
use std::collections::HashMap;

use crate::encoding::Symbol;
use crate::references::{HaplotypeRef, HaplotypeRefTrait};

#[cfg(not(feature = "parallel"))]
pub type Store<S> = HashMap<usize, HaplotypeRef<S>>;
#[cfg(feature = "parallel")]
pub type Store<S> = DashMap<usize, HaplotypeRef<S>>;

pub trait HaplotypeStore:
    Clone + std::fmt::Debug + Default + FromIterator<(usize, Self::Item)>
{
    type Item: HaplotypeRefTrait;
    type Symbol: Symbol;

    fn get_value(&self, id: &usize) -> Option<Self::Item>;
    fn put_value(&mut self, id: usize, haplotype: Self::Item);
    fn extend_store(&mut self, haplotypes: Self);
}

pub trait SyncHaplotypeStore: HaplotypeStore + Send + Sync {
    type Item: HaplotypeRefTrait;
    type Symbol: Symbol;

    fn put_value_sync(&self, id: usize, haplotype: <Self as SyncHaplotypeStore>::Item);
}

impl<S: Symbol> HaplotypeStore for HashMap<usize, HaplotypeRef<S>> {
    type Item = HaplotypeRef<S>;
    type Symbol = S;

    fn get_value(&self, id: &usize) -> Option<Self::Item> {
        self.get(id).cloned()
    }

    fn put_value(&mut self, id: usize, haplotype: Self::Item) {
        self.entry(id).or_insert(haplotype);
    }

    fn extend_store(&mut self, haplotypes: Self) {
        std::iter::Extend::extend(self, haplotypes);
    }
}

#[cfg(feature = "parallel")]
impl<S: Symbol> HaplotypeStore for DashMap<usize, HaplotypeRef<S>> {
    type Item = HaplotypeRef<S>;
    type Symbol = S;

    fn get_value(&self, id: &usize) -> Option<HaplotypeRef<S>> {
        self.get(id).map(|entry| entry.clone())
    }

    fn put_value(&mut self, id: usize, haplotype: HaplotypeRef<S>) {
        self.entry(id).or_insert(haplotype);
    }

    fn extend_store(&mut self, haplotypes: Self) {
        for (id, haplotype) in haplotypes {
            self.insert(id, haplotype);
        }
    }
}

#[cfg(feature = "parallel")]
impl<S: Symbol> SyncHaplotypeStore for DashMap<usize, HaplotypeRef<S>> {
    type Item = HaplotypeRef<S>;
    type Symbol = S;

    fn put_value_sync(&self, id: usize, haplotype: HaplotypeRef<S>) {
        self.entry(id).or_insert(haplotype);
    }
}

/// Creates a `Population` from the arguments.
///
/// `population!` creates a `Population` with the same syntax as array expressions.
///
/// - Create a `Population` with a list of haplotypes:
///
/// ```rust
/// # use virolution::encoding::Nucleotide as Nt;
/// # use virolution::core::AttributeSetDefinition;
/// # use virolution::core::haplotype::Wildtype;
/// # use virolution::core::population::{Store, Population};
/// # use virolution::population;
/// let wt = Wildtype::new(vec![Nt::A; 10], &AttributeSetDefinition::new());
/// let ht = wt.create_descendant(vec![0], vec![Nt::C]);
/// let population: Population<Store<Nt>> = population![&wt; &ht];
/// ```
///
/// - Create a `Population` with a haplotype and a size:
///
/// ```rust
/// # use virolution::encoding::Nucleotide as Nt;
/// # use virolution::core::AttributeSetDefinition;
/// # use virolution::core::haplotype::Wildtype;
/// # use virolution::core::population::{Store, Population};
/// # use virolution::population;
/// let wt = Wildtype::new(vec![Nt::A; 10], &AttributeSetDefinition::new());
/// let population: Population<Store<Nt>> = population![wt, 10];
/// ```
///
/// - Create a `Population` from a list of haplotypes and sizes:
///
/// ```rust
/// # use virolution::encoding::Nucleotide as Nt;
/// # use virolution::core::AttributeSetDefinition;
/// # use virolution::core::haplotype::Wildtype;
/// # use virolution::core::population::{Store, Population};
/// # use virolution::population;
/// let wt = Wildtype::new(vec![Nt::A; 10], &AttributeSetDefinition::new());
/// let ht = wt.create_descendant(vec![0], vec![Nt::C]);
/// let population: Population<Store<Nt>> = population![wt, 10; ht, 5];
/// ```
///
#[macro_export]
macro_rules! population {
    () => {
        $crate::core::Population::new()
    };
    ($haplotype:expr, $size:expr) => {
        $crate::core::Population::from_haplotype($haplotype, $size)
    };
    ($( $haplotype:expr );+) => {
        {
        let mut population = $crate::core::Population::new();
        $(
            population.push($haplotype);
        )+
        population
    }
    };
    ($( $haplotype:expr , $size:expr );+) => {
        $crate::core::Population::from_iter(vec![$( $crate::core::Population::from_haplotype($haplotype, $size) ),+])
    };
}

/// A `Population` is a collection of haplotypes.
#[derive(Clone, Default, Debug)]
pub struct Population<M: HaplotypeStore> {
    population: Vec<Cell<usize>>,
    haplotypes: M,
    // we need to store the generic type S to use it in haplotype store. PhantomData is used to
    // store the type during compile time, but it is not used during runtime.
    // _marker: std::marker::PhantomData<S>,
}

unsafe impl<M: HaplotypeStore + Send> Send for Population<M> {}
unsafe impl<M: HaplotypeStore + Sync> Sync for Population<M> {}

pub struct PopulationIterator<'a, M: HaplotypeStore> {
    population: &'a Population<M>,
    index: usize,
}

impl<M: HaplotypeStore> Iterator for PopulationIterator<'_, M> {
    type Item = M::Item;

    fn next(&mut self) -> Option<M::Item> {
        if self.index >= self.population.len() {
            return None;
        }

        let ref_id = self.population.population[self.index].get();
        self.index += 1;
        self.population
            .haplotypes
            .get_value(&ref_id)
            .or_else(|| panic!("No haplotype with index {ref_id}"))
    }
}

impl<M: HaplotypeStore> FromIterator<Population<M>> for Population<M> {
    fn from_iter<I: IntoIterator<Item = Self>>(iter: I) -> Self {
        let mut population: Vec<Cell<usize>> = Vec::new();
        let mut haplotypes: M = M::default();

        for pop in iter {
            population.extend(pop.population);
            haplotypes.extend_store(pop.haplotypes);
        }

        Self {
            population,
            haplotypes,
        }
    }
}

impl<M: HaplotypeStore> PartialEq for Population<M> {
    fn eq(&self, other: &Self) -> bool {
        self.population == other.population
    }
}

#[cfg(feature = "parallel")]
impl<M> Population<M>
where
    M: SyncHaplotypeStore,
{
    /// Insert a `HaplotypeRef` at a specific position.
    ///
    /// This method is partially thread safe, if it it is operated on disjoint positions of the
    /// population. If the same position is accessed by multiple threads, the behavior is
    /// undefined.
    pub fn insert_sync(&self, position: &usize, haplotype: <M as SyncHaplotypeStore>::Item) {
        let ref_id = haplotype.get_id();
        self.population[*position].set(ref_id);
        self.haplotypes.put_value_sync(ref_id, haplotype);
    }

    /// Get a random subsample `Population` of the `Population` with parallel sampling.
    pub fn par_sample(&self, size: usize, weights: &[usize]) -> Self {
        let sampler = WeightedAliasIndex::new(weights.to_vec()).unwrap();

        let population: Vec<Cell<usize>> = (0..size)
            .into_par_iter()
            .map_init(rand::rng, |rng, _| {
                self.population[sampler.sample(rng)].clone()
            })
            .collect();
        let haplotypes: M = population
            .iter()
            .map(|cell| cell.get())
            .unique()
            .map(|id| {
                (
                    id,
                    self.haplotypes
                        .get_value(&id)
                        .expect("Haplotype not found in store.")
                        .clone(),
                )
            })
            .collect();

        Self {
            population,
            haplotypes,
        }
    }
}

impl<M: HaplotypeStore> Population<M> {
    /// Construct a new, empty `Population`.
    pub fn new() -> Self {
        Self {
            population: Vec::new(),
            haplotypes: Default::default(),
        }
    }

    /// Construct a new `Population` from a `HaplotypeRef` and a size.
    pub fn from_haplotype(haplotype: M::Item, size: usize) -> Self {
        let ref_id = haplotype.get_id();
        let population = (0..size).map(|_| Cell::new(ref_id)).collect();
        let mut haplotypes = M::default();
        haplotypes.put_value(haplotype.get_id(), haplotype);
        Self {
            population,
            haplotypes,
        }
    }

    /// Construct a new `Population` from a `Vec` of `HaplotypeRef`s.
    pub fn from_references(references: Vec<M::Item>) -> Self {
        let mut population = Vec::new();
        let mut haplotypes = M::default();

        for reference in references {
            let ref_id = reference.get_id();
            population.push(Cell::new(ref_id));
            haplotypes.put_value(ref_id, reference);
        }

        Self {
            population,
            haplotypes,
        }
    }

    /// Get an iterator over the `Population`.
    pub fn iter(&self) -> PopulationIterator<'_, M> {
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
    pub fn insert(&mut self, position: &usize, haplotype: M::Item) {
        let ref_id = haplotype.get_id();
        self.population[*position] = ref_id.into();
        self.haplotypes.put_value(ref_id, haplotype);
    }

    /// Push a `HaplotypeRef` to the end of the `Population`.
    pub fn push(&mut self, haplotype: &M::Item) {
        let ref_id = haplotype.get_id();
        self.population.push(ref_id.into());
        self.haplotypes.put_value(ref_id, haplotype.clone());
    }

    /// Get a reference to the `Population` vector.
    pub fn get_population(&self) -> &Vec<Cell<usize>> {
        &self.population
    }

    /// Get a reference to the `Haplotypes` map.
    pub fn get_haplotypes(&self) -> &M {
        &self.haplotypes
    }

    /// Get a `HaplotypeRef` from the `Population` by index.
    pub fn get(&self, index: &usize) -> M::Item {
        let ref_id = self.population[*index].get();
        self.haplotypes
            .get_value(&ref_id)
            .unwrap_or_else(|| panic!("No haplotype with index {index}"))
    }

    /// Choose `Population` of size `size` from the `Population`.
    pub fn choose_multiple(&self, rng: &mut ThreadRng, size: usize) -> Self {
        let population: Vec<Cell<usize>> = self
            .population
            .choose_multiple(rng, size)
            .cloned()
            .collect();
        let haplotypes = population
            .iter()
            .map(|cell| cell.get())
            .unique()
            .map(|id| (id, self.haplotypes.get_value(&id).unwrap().clone()))
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
            .map(|cell| cell.get())
            .unique()
            .map(|id| (id, self.haplotypes.get_value(&id).unwrap().clone()))
            .collect();
    }

    /// Get a random subsample `Population` of the `Population`.
    pub fn sample(&self, size: usize, weights: &[usize]) -> Self {
        let sampler = WeightedAliasIndex::new(weights.to_vec()).unwrap();

        let mut rng = rand::rng();
        let population: Vec<Cell<usize>> = (0..size)
            .map(|_| self.population[sampler.sample(&mut rng)].clone())
            .collect();
        let haplotypes: M = population
            .iter()
            .map(|cell| cell.get())
            .unique()
            .map(|id| (id, self.haplotypes.get_value(&id).unwrap().clone()))
            .collect();

        Self {
            population,
            haplotypes,
        }
    }

    /// Select a `Population` from the `Population` based on a vector of indices.
    pub fn select(&self, indices: &[usize]) -> Self {
        let population: Vec<Cell<usize>> = indices
            .iter()
            .map(|&index| self.population[index].clone())
            .collect();
        let haplotypes: M = population
            .iter()
            .map(|cell| cell.get())
            .unique()
            .map(|id| {
                (
                    id,
                    self.haplotypes
                        .get_value(&id)
                        .expect("Haplotype not found in store.")
                        .clone(),
                )
            })
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

    use crate::core::attributes::AttributeSetDefinition;
    use crate::core::haplotype::Wildtype;
    use crate::encoding::Nucleotide as Nt;

    #[test]
    fn is_empty() {
        let mut population: Population<Store<Nt>> = Population::new();
        assert!(population.is_empty());
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt);
        assert!(!population.is_empty());
    }

    #[test]
    fn len() {
        let mut population: Population<Store<Nt>> = Population::new();
        assert_eq!(population.len(), 0);

        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt);
        assert_eq!(population.len(), 1);
    }

    #[test]
    fn insert() {
        let mut population: Population<Store<Nt>> = Population::new();
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt);
        assert_eq!(population.len(), 1);
        assert_eq!(population.get(&0), wt);

        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.insert(&0, wt2.clone());
        assert_eq!(population.len(), 1);
        assert_ne!(population.get(&0), wt);
        assert_eq!(population.get(&0), wt2);
    }

    #[test]
    fn push() {
        let mut population: Population<Store<Nt>> = Population::new();
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt);
        assert_eq!(population.len(), 1);
        assert_eq!(population.get(&0), wt);

        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt2);
        assert_eq!(population.len(), 2);
        assert_eq!(population.get(&1), wt2);
    }

    #[test]
    fn iterate() {
        let mut population: Population<Store<Nt>> = Population::new();
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt);
        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        population.push(&wt2);

        let mut iter = population.iter();
        assert_eq!(iter.next().unwrap().get_id(), wt.get_id());
        assert_eq!(iter.next().unwrap().get_id(), wt2.get_id());
        assert_eq!(iter.next(), None);
    }

    #[test]
    fn from_iter() {
        let attribute_definition = AttributeSetDefinition::new();
        let wt1 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let mut population1: Population<Store<Nt>> = Population::new();
        population1.push(&wt1);

        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let mut population2 = Population::new();
        population2.push(&wt2);

        let population = Population::from_iter(vec![population1, population2]);
        assert_eq!(population.len(), 2);
        assert_eq!(population.get(&0), wt1);
        assert_eq!(population.get(&1), wt2);
    }

    #[test]
    fn macro_empty() {
        let population: Population<Store<Nt>> = population![];
        assert_eq!(population.len(), 0);
        assert!(population.is_empty());
    }

    #[test]
    fn macro_from_haplotype() {
        let attribute_definition = AttributeSetDefinition::new();
        let wt1 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let pop1: Population<Store<Nt>> = population![wt1.clone(), 1];
        assert_eq!(pop1.len(), 1);
        let pop2: Population<Store<Nt>> = population![wt1.clone(), 2];
        assert_eq!(pop2.len(), 2);
        let pop3: Population<Store<Nt>> = population![wt1.clone(), 10];
        assert_eq!(pop3.len(), 10);
        let pop4: Population<Store<Nt>> = population![wt1.clone(), 1_000_000];
        assert_eq!(pop4.len(), 1_000_000);
    }

    #[test]
    fn macro_from_haplotypes() {
        let attribute_definition = AttributeSetDefinition::new();
        let wt1 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let wt3 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let population: Population<Store<Nt>> = population![&wt1; &wt2; &wt3];

        assert_eq!(population.len(), 3);
        assert_eq!(population.get(&0), wt1);
        assert_eq!(population.get(&1), wt2);
        assert_eq!(population.get(&2), wt3);
    }

    #[test]
    fn macro_from_haplotypes_variable_lengths() {
        let attribute_definition = AttributeSetDefinition::new();
        let wt1 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let wt2 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let wt3 = Wildtype::new(vec![Nt::A; 10], &attribute_definition);
        let population: Population<Store<Nt>> =
            population![wt1.clone(), 2; wt2.clone(), 1; wt3.clone(), 3];

        assert_eq!(population.len(), 6);
        assert_eq!(population.get(&0), wt1);
        assert_eq!(population.get(&1), wt1);
        assert_eq!(population.get(&2), wt2);
        assert_eq!(population.get(&3), wt3);
        assert_eq!(population.get(&4), wt3);
        assert_eq!(population.get(&5), wt3);
    }
}
