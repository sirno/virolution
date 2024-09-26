//! Haplotype representation and operations
//!
//! The `Haplotype` struct represents a single haplotype, represented as part of
//! a tree structure. The tree structure is implemented using weak references
//! and strong backreferences allowing any resources that are not observable
//! anymore to be dropped.
//!
//! The `Haplotype` struct is an enum that can be a representation of any
//! possible transformation of a preceeding haplotype. Currently, the following
//! transformations are supported:
//!
//! - Wildtype: The wildtype haplotype is the starting point of any phylogeny.
//!   It is the only haplotype that does not have a parent.
//! - Mutant: A mutant haplotype is a haplotype that has a set of mutations that
//!   are not present in its parent. The mutations are stored as a HashMap of
//!   positions and symbols.
//! - Recombinant: A recombinant haplotype is a haplotype that consists of a
//!   combination of two parent haplotypes.
//!

use derivative::Derivative;
use seq_io::fasta::OwnedRecord;
use std::collections::HashMap;
use std::fmt;
use std::sync::atomic::{AtomicIsize, Ordering};
use std::sync::OnceLock;

use crate::encoding::STRICT_ENCODE;
use crate::references::DescendantsCell;
use crate::references::{HaplotypeRef, HaplotypeWeak};

use super::cache::{CachedValue, VirolutionCache};
use super::fitness::FitnessProvider;

pub static N_FITNESS_TABLES: OnceLock<usize> = OnceLock::new();

pub fn set_number_of_fitness_tables(value: usize) -> Result<(), usize> {
    if N_FITNESS_TABLES.get().is_some() {
        eprintln!("Trying to set `N_FITNESS_TABLES`, which has already been set.");
    }
    N_FITNESS_TABLES.set(value)
}

fn make_fitness_cache() -> Vec<OnceLock<f64>> {
    // warn if the number of fitness providers has not been explicitly set at runtime
    if N_FITNESS_TABLES.get().is_none() {
        eprintln!("`N_FITNESS_TABLES` is not set, using default value of 1.");
    }
    std::iter::repeat_with(OnceLock::new)
        .take(*N_FITNESS_TABLES.get_or_init(|| 1))
        .collect()
}

// #[derive(Clone, Debug, Deref)]
pub type Symbol = Option<u8>;

#[derive(Debug)]
pub enum Haplotype {
    Wildtype(Wildtype),
    Mutant(Mutant),
    Recombinant(Recombinant),
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Wildtype {
    reference: HaplotypeWeak,
    sequence: Vec<Symbol>,
    descendants: DescendantsCell,
    dirty_descendants: AtomicIsize,
}

#[derive(Debug)]
pub struct Mutant {
    reference: HaplotypeWeak,
    wildtype: HaplotypeWeak,
    ancestor: HaplotypeRef,
    changes: HashMap<usize, (Symbol, Symbol)>,
    generation: usize,
    fitness: Vec<OnceLock<f64>>,
    descendants: DescendantsCell,
    dirty_descendants: AtomicIsize,
}

#[derive(Debug)]
pub struct Recombinant {
    reference: HaplotypeWeak,
    wildtype: HaplotypeWeak,
    left_ancestor: HaplotypeRef,
    right_ancestor: HaplotypeRef,
    left_position: usize,
    right_position: usize,
    generation: usize,
    fitness: Vec<OnceLock<f64>>,
    descendants: DescendantsCell,
    dirty_descendants: AtomicIsize,
}

// We have to explicitly implement Send and Sync for Haplotype, because it contains several
// read-only fields that are not thread-safe. For example, it may contain a HashMap, which reflects
// the changes in the haplotype data, that is created during the construction of the Haplotype.
// To my knowledge the only effect of these lines (right now) is that it silences the warning
// that Arc references an object that is not Send and Sync.
#[cfg(feature = "parallel")]
unsafe impl Send for Haplotype {}
#[cfg(feature = "parallel")]
unsafe impl Sync for Haplotype {}

impl Drop for Haplotype {
    fn drop(&mut self) {
        match self {
            Haplotype::Wildtype(_wt) => {}
            Haplotype::Mutant(mt) => {
                mt.ancestor.increment_dirty_descendants();
            }
            Haplotype::Recombinant(rc) => {
                rc.left_ancestor.increment_dirty_descendants();
                rc.right_ancestor.increment_dirty_descendants();
            }
        }
    }
}

impl Haplotype {
    pub fn create_descendant(
        &self,
        positions: Vec<usize>,
        changes: Vec<Symbol>,
        generation: usize,
    ) -> HaplotypeRef {
        let ancestor = self.get_reference();
        let wildtype = self.get_wildtype();

        let changes = HashMap::from_iter(
            positions
                .iter()
                .zip(changes.iter())
                .map(|(pos, sym)| (*pos, (self.get_base(pos), *sym))),
        );

        changes.iter().for_each(|(pos, change)| {
            if change.0 == change.1 {
                dbg!(pos, change);
            }
        });

        let descendant = Mutant::new(ancestor, wildtype.get_weak(), changes, generation);

        self.add_descendant(descendant.get_weak());

        descendant
    }

    pub fn create_recombinant(
        left_ancestor: &HaplotypeRef,
        right_ancestor: &HaplotypeRef,
        left_position: usize,
        right_position: usize,
        generation: usize,
    ) -> HaplotypeRef {
        let wildtype = left_ancestor.get_wildtype();

        let recombinant = Recombinant::new(
            wildtype,
            left_ancestor.clone(),
            right_ancestor.clone(),
            left_position,
            right_position,
            generation,
        );

        left_ancestor.add_descendant(recombinant.get_weak());
        right_ancestor.add_descendant(recombinant.get_weak());
        recombinant
    }

    pub(crate) fn is_mutant(&self) -> bool {
        matches!(self, Haplotype::Mutant(_))
    }

    /// Unwraps the haplotype into a mutant.
    ///
    /// This function will panic if the haplotype is not a mutant and is only intended for internal
    /// use while minimizing the tree.
    pub(crate) fn try_unwrap_mutant(&self) -> Option<&Mutant> {
        match self {
            Haplotype::Mutant(ht) => Some(ht),
            _ => None,
        }
    }

    /// Returns a reference to the changes that are present in the haplotype if the type allows it.
    pub fn try_get_changes(&self) -> Option<&HashMap<usize, (Symbol, Symbol)>> {
        match self {
            Haplotype::Mutant(ht) => Some(&ht.changes),
            _ => None,
        }
    }

    /// Returns a reference to an ancestor if the type allows it.
    ///
    /// If there are multiple ancestors, in the future any ancestor may be returned.
    pub fn try_get_ancestor(&self) -> Option<&HaplotypeRef> {
        match self {
            Haplotype::Mutant(ht) => Some(&ht.ancestor),
            Haplotype::Recombinant(rc) => Some(&rc.left_ancestor),
            _ => None,
        }
    }

    pub fn get_reference(&self) -> HaplotypeRef {
        let weak = match self {
            Haplotype::Wildtype(wt) => &wt.reference,
            Haplotype::Mutant(ht) => &ht.reference,
            Haplotype::Recombinant(rc) => &rc.reference,
        };
        weak.upgrade().expect("Self-reference has been dropped.")
    }

    pub fn get_wildtype(&self) -> HaplotypeRef {
        match self {
            Haplotype::Wildtype(wt) => wt.get_reference().clone(),
            Haplotype::Mutant(ht) => ht.wildtype.upgrade().unwrap().clone(),
            Haplotype::Recombinant(rc) => rc.wildtype.upgrade().unwrap().clone(),
        }
    }

    pub fn get_ancestors(&self) -> (Option<HaplotypeRef>, Option<HaplotypeRef>) {
        match self {
            Haplotype::Mutant(ht) => (Some(ht.ancestor.clone()), None),
            Haplotype::Recombinant(rc) => (
                Some(rc.left_ancestor.clone()),
                Some(rc.right_ancestor.clone()),
            ),
            _ => (None, None),
        }
    }

    pub fn get_descendants(&self) -> &DescendantsCell {
        match self {
            Haplotype::Wildtype(wt) => &wt.descendants,
            Haplotype::Mutant(ht) => &ht.descendants,
            Haplotype::Recombinant(rc) => &rc.descendants,
        }
    }

    /// Returns the number of descendants.
    fn add_descendant(&self, descendant: HaplotypeWeak) {
        let (descendants, dirty_descendants) = match self {
            Haplotype::Wildtype(wt) => (&wt.descendants, &wt.dirty_descendants),
            Haplotype::Mutant(ht) => (&ht.descendants, &ht.dirty_descendants),
            Haplotype::Recombinant(rc) => (&rc.descendants, &rc.dirty_descendants),
        };

        let mut descendants_guard = descendants.lock();

        // if there are dirty descendants, replace one of them
        if dirty_descendants.load(Ordering::Relaxed) > 0
            && let Some(idx) = descendants_guard.iter().position(|x| !x.exists())
        {
            descendants_guard[idx] = descendant;
            dirty_descendants.fetch_add(1, Ordering::Relaxed);
            return;
        }

        descendants_guard.push(descendant);
    }

    pub(crate) fn increment_dirty_descendants(&self) {
        let dirty_descendants = match self {
            Haplotype::Wildtype(wt) => &wt.dirty_descendants,
            Haplotype::Mutant(ht) => &ht.dirty_descendants,
            Haplotype::Recombinant(rc) => &rc.dirty_descendants,
        };
        dirty_descendants.fetch_sub(1, Ordering::Relaxed);
    }

    pub fn get_base(&self, position: &usize) -> Symbol {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Mutant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_generation(&self) -> usize {
        match self {
            Haplotype::Wildtype(_wt) => 0,
            Haplotype::Mutant(ht) => ht.get_generation(),
            Haplotype::Recombinant(rc) => rc.get_generation(),
        }
    }

    pub fn get_sequence(&self) -> Vec<Symbol> {
        let mutations = self.get_mutations();
        let mut sequence = self.get_wildtype_sequence();

        for (&position, &new) in mutations.iter() {
            sequence[position] = new;
        }

        sequence
    }

    pub fn get_wildtype_sequence(&self) -> Vec<Symbol> {
        match self {
            Haplotype::Wildtype(wt) => wt.sequence.to_vec(),
            Haplotype::Mutant(_ht) => self.get_wildtype().get_wildtype_sequence(),
            Haplotype::Recombinant(_rc) => self.get_wildtype().get_wildtype_sequence(),
        }
    }

    pub fn get_string(&self) -> String {
        let mutations = self.get_mutations();
        let wildtype = self.get_wildtype_sequence();

        let mut out = String::new();
        for (position, to) in mutations.iter() {
            let from = wildtype[*position];
            match (from, to) {
                (Some(f), Some(t)) => {
                    out.push_str(
                        format!(
                            ";{position}:{}->{}",
                            char::from(STRICT_ENCODE[&f]),
                            char::from(STRICT_ENCODE[t])
                        )
                        .as_str(),
                    );
                }
                (None, Some(t)) => {
                    if let Some(wt_symbol) = wildtype[*position] {
                        out.push_str(
                            format!(
                                ";{position}:{}->{}",
                                char::from(STRICT_ENCODE[&wt_symbol]),
                                char::from(STRICT_ENCODE[t])
                            )
                            .as_str(),
                        )
                    }
                }
                _ => {}
            }
        }

        if out.is_empty() {
            return "wt".to_string();
        }

        out.remove(0);
        out
    }

    /// Returns a HashMap of mutations that are present in the haplotype.
    ///
    /// Calling this with large trees can be expensive, especially when there are many
    /// recombinants. For every recombinant, both ancestors have to be traversed to find
    /// the mutations that are present in the recombinant.
    pub fn get_mutations(&self) -> CachedValue<HashMap<usize, Symbol>> {
        match self {
            Haplotype::Wildtype(_wt) => Default::default(),
            Haplotype::Mutant(ht) => ht.get_mutations(),
            Haplotype::Recombinant(rc) => rc.get_mutations(),
        }
    }

    pub fn get_raw_fitness(&self, fitness_provider: &FitnessProvider) -> f64 {
        match self {
            Haplotype::Wildtype(_wt) => 1.,
            Haplotype::Mutant(ht) => ht.get_raw_fitness(fitness_provider),
            Haplotype::Recombinant(rc) => rc.get_raw_fitness(fitness_provider),
        }
    }

    pub fn get_fitness(&self, fitness_provider: &FitnessProvider) -> f64 {
        match self {
            Haplotype::Wildtype(_wt) => 1.,
            Haplotype::Mutant(ht) => ht.get_fitness(fitness_provider),
            Haplotype::Recombinant(rc) => rc.get_fitness(fitness_provider),
        }
    }

    pub fn get_record(&self, head: &str) -> OwnedRecord {
        OwnedRecord {
            head: head.to_string().as_bytes().to_vec(),
            seq: self
                .get_sequence()
                .into_iter()
                .map(|symbol| match symbol {
                    Some(s) => STRICT_ENCODE[&s],
                    None => 0x2d,
                })
                .collect(),
        }
    }

    pub fn get_length(&self) -> usize {
        match self {
            Haplotype::Wildtype(wt) => wt.get_length(),
            Haplotype::Mutant(ht) => ht.get_length(),
            Haplotype::Recombinant(rc) => rc.get_length(),
        }
    }

    pub fn get_tree(&self) -> String {
        let tree = self.get_subtree(self.get_reference().get_weak());
        format!("{};", tree)
    }

    pub fn get_subtree(&self, ancestor: HaplotypeWeak) -> String {
        match self {
            Haplotype::Wildtype(wt) => wt.get_subtree(),
            Haplotype::Mutant(ht) => ht.get_subtree(),
            Haplotype::Recombinant(rc) => rc.get_subtree(ancestor),
        }
    }
}

impl Wildtype {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(sequence: Vec<Symbol>) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Wildtype(Self {
                reference: reference.clone(),
                sequence: sequence.clone(),
                descendants: DescendantsCell::new(),
                dirty_descendants: AtomicIsize::new(0),
            })
        })
    }

    pub fn get_reference(&self) -> HaplotypeRef {
        self.reference
            .upgrade()
            .expect("Self-reference has been dropped.")
    }

    pub fn get_base(&self, position: &usize) -> Symbol {
        self.sequence[*position]
    }

    pub fn get_length(&self) -> usize {
        self.sequence.len()
    }

    pub fn get_subtree(&self) -> String {
        let inner = self
            .descendants
            .lock()
            .iter()
            .filter_map(|x| x.upgrade())
            .map(|x| x.get_subtree(self.reference.clone()))
            .collect::<Vec<String>>()
            .join(",");
        format!("({inner})wt")
    }
}

impl Mutant {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        ancestor: HaplotypeRef,
        wildtype: HaplotypeWeak,
        changes: HashMap<usize, (Symbol, Symbol)>,
        generation: usize,
    ) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Mutant(Self {
                reference: reference.clone(),
                wildtype: wildtype.clone(),
                ancestor: ancestor.clone(),
                changes: changes.clone(),
                generation,
                fitness: make_fitness_cache(),
                descendants: DescendantsCell::new(),
                dirty_descendants: AtomicIsize::new(0),
            })
        })
    }

    /// Creates a new mutant haplotype and replaces `last` with it.
    ///
    /// This function is not thread-safe and should only be called by a single thread at a time.
    pub(crate) unsafe fn new_and_replace(
        last: &Mutant,
        ancestor: HaplotypeRef,
        wildtype: HaplotypeWeak,
        changes: HashMap<usize, (Symbol, Symbol)>,
        generation: usize,
    ) {
        // collect all descendants that are still alive
        let descendants: Vec<HaplotypeWeak> = last
            .descendants
            .lock()
            .iter()
            .filter(|x| x.exists())
            .cloned()
            .collect();

        // make the ancestor aware of the new descendant
        ancestor.add_descendant(last.reference.clone());

        // create new node
        let new = HaplotypeRef::new(Haplotype::Mutant(Self {
            reference: last.reference.clone(),
            wildtype,
            ancestor,
            changes,
            generation,
            fitness: make_fitness_cache(),
            descendants: DescendantsCell::from_iter(descendants),
            dirty_descendants: AtomicIsize::new(0),
        }));

        let old_ptr = last.reference.as_ptr() as *mut Haplotype;
        let new_ptr = new.as_ptr() as *mut Haplotype;

        // replace the old reference with the new one
        // this swaps the reference count of the old reference for the new one
        std::ptr::swap(old_ptr, new_ptr);
    }

    pub fn get_base(&self, position: &usize) -> Symbol {
        match self.changes.get(position) {
            Some((_from, to)) => *to,
            None => self.ancestor.get_base(position),
        }
    }

    pub fn iter_changes(&self) -> impl Iterator<Item = (&usize, &(Symbol, Symbol))> + '_ {
        self.changes.iter()
    }

    pub fn get_generation(&self) -> usize {
        self.generation
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    pub fn get_ancestor(&self) -> HaplotypeRef {
        self.ancestor.clone()
    }

    pub fn get_changes(&self) -> &HashMap<usize, (Symbol, Symbol)> {
        &self.changes
    }

    pub fn get_mutations(&self) -> CachedValue<HashMap<usize, Symbol>> {
        let mut mutations = self.ancestor.get_mutations().clone_inner();

        let wt_ref = self.wildtype.upgrade().unwrap();

        self.changes.iter().for_each(|(position, (_, new))| {
            let wt_base = wt_ref.get_base(position);
            if *new == wt_base {
                mutations.remove(position);
            } else {
                mutations.insert(*position, *new);
            }
        });

        CachedValue::new(mutations)
    }

    pub fn get_raw_fitness(&self, fitness_provider: &FitnessProvider) -> f64 {
        *self
            .fitness
            .get(fitness_provider.id)
            .unwrap()
            .get_or_init(|| fitness_provider.get_fitness(&self.reference.upgrade().unwrap()))
    }

    pub fn get_fitness(&self, fitness_provider: &FitnessProvider) -> f64 {
        fitness_provider.apply_utility(self.get_raw_fitness(fitness_provider))
    }

    pub fn get_subtree(&self) -> String {
        let current = self.reference.upgrade().unwrap();

        let block_id = current.get_block_id();
        let branch_length = current.get_generation() - self.ancestor.get_generation();

        let node = format!("'{block_id}':{branch_length}");

        let descendants = current.get_descendants();
        let descendants_guard = descendants.lock();

        let descendants = descendants_guard
            .iter()
            .filter_map(|x| x.upgrade())
            .map(|x| x.get_subtree(current.get_weak()))
            .collect::<Vec<String>>()
            .join(",");

        if descendants.is_empty() {
            return node;
        }

        format!("({descendants}){node}")
    }

    pub(crate) fn try_merge_node(&self) {
        let n_descendants = self
            .descendants
            .lock()
            .len()
            .checked_add_signed(self.dirty_descendants.load(Ordering::Relaxed))
            .expect("Number of descendants overflowed.");

        let descendant = match n_descendants {
            1 => self
                .descendants
                .lock()
                .iter()
                .filter_map(|x| x.upgrade())
                .next(),
            _ => None,
        };

        let merger: [Option<&Mutant>; 3] = [
            descendant.as_ref().and_then(|x| x.try_unwrap_mutant()),
            Some(self),
            self.ancestor.try_unwrap_mutant(),
        ];

        // if there is nothing to merge, return
        if merger[0].is_none() && merger[2].is_none() {
            return;
        }

        // aggregate changes
        let changes: HashMap<usize, (Symbol, Symbol)> = merger
            .iter()
            .rev()
            .flatten()
            .flat_map(|x| x.changes.iter())
            .map(|(position, change)| (*position, *change))
            .collect();

        // last is the node furthest away from the wildtype
        // this node should be replaced, it is the node that may still be
        // visible
        let last = *merger.iter().flatten().next().unwrap();

        // first is the node closest to the wildtype
        let first = *merger.iter().rev().flatten().next().unwrap();

        // determine generation
        let generation = last.get_generation();

        // extract ancestor and wildtype
        let ancestor = first.ancestor.clone();
        let wildtype = first.wildtype.clone();

        // // we need to decrement the dirty descendants of the ancestor
        // // because it will be decremented when the previous node gets removed
        // ancestor.decrement_dirty_descendants();

        // create new node
        unsafe {
            Mutant::new_and_replace(last, ancestor, wildtype, changes, generation);
        };
    }
}

impl Recombinant {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        wildtype: HaplotypeRef,
        left_ancestor: HaplotypeRef,
        right_ancestor: HaplotypeRef,
        left_position: usize,
        right_position: usize,
        generation: usize,
    ) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Recombinant(Self {
                reference: reference.clone(),
                wildtype: wildtype.get_weak(),
                left_ancestor: left_ancestor.clone(),
                right_ancestor: right_ancestor.clone(),
                left_position,
                right_position,
                generation,
                fitness: make_fitness_cache(),
                descendants: DescendantsCell::new(),
                dirty_descendants: AtomicIsize::new(0),
            })
        })
    }

    pub fn get_base(&self, position: &usize) -> Symbol {
        if *position >= self.left_position && *position < self.right_position {
            return self.left_ancestor.get_base(position);
        }

        self.right_ancestor.get_base(position)
    }

    pub fn get_generation(&self) -> usize {
        self.generation
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    fn get_mutations_cache(
    ) -> &'static ::cached::once_cell::sync::Lazy<VirolutionCache<HashMap<usize, Symbol>>> {
        static MUTATIONS_CACHE: ::cached::once_cell::sync::Lazy<
            VirolutionCache<HashMap<usize, Symbol>>,
        > = ::cached::once_cell::sync::Lazy::new(|| VirolutionCache::new(100));
        &MUTATIONS_CACHE
    }

    pub fn get_mutations(&self) -> CachedValue<HashMap<usize, Symbol>> {
        let key = self.reference.get_id();

        // try to use the cache first
        let cache = Self::get_mutations_cache();
        if let Some(mutations) = cache.cache_get(&key) {
            return mutations;
        }

        // collect all mutations
        let left_mutations = self.left_ancestor.get_mutations();
        let right_mutations = self.right_ancestor.get_mutations();

        let mut mutations = HashMap::new();

        for (&position, &change) in left_mutations.iter() {
            if position >= self.left_position && position < self.right_position {
                mutations.insert(position, change);
            }
        }

        for (&position, &change) in right_mutations.iter() {
            if position < self.left_position || position >= self.right_position {
                mutations.insert(position, change);
            }
        }

        // set and return computed value
        cache.cache_set(key, mutations.clone())
    }

    pub fn get_raw_fitness(&self, fitness_provider: &FitnessProvider) -> f64 {
        *self
            .fitness
            .get(fitness_provider.id)
            .unwrap()
            .get_or_init(|| fitness_provider.get_fitness(&self.reference.upgrade().unwrap()))
    }

    pub fn get_fitness(&self, fitness_provider: &FitnessProvider) -> f64 {
        fitness_provider.apply_utility(self.get_raw_fitness(fitness_provider))
    }

    pub fn get_subtree(&self, ancestor: HaplotypeWeak) -> String {
        let block_id = self.reference.get_block_id();
        let branch_length = self.generation - ancestor.upgrade().unwrap().get_generation();

        let node = format!("#R'{block_id}':{branch_length}");

        let descendants = if self.left_ancestor.get_weak() == ancestor {
            self.descendants
                .lock()
                .iter()
                .filter_map(|x| x.upgrade())
                .map(|x| x.get_subtree(self.reference.clone()))
                .collect::<Vec<String>>()
                .join(",")
        } else {
            "".to_string()
        };

        if descendants.is_empty() {
            return node;
        }

        format!("({descendants}){node}",)
    }
}

impl fmt::Display for Haplotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.get_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serial_test::serial;

    #[test]
    fn initiate_wildtype() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);
        for i in 0..4 {
            assert_eq!(wt.get_base(&i), Some(i as u8));
        }
    }

    #[test]
    fn create_descendant() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);
        let ht = wt.create_descendant(vec![0], vec![Some(0x03)], 0);
        assert_eq!(ht.get_base(&0), Some(0x03));
        assert_eq!(ht.get_base(&1), Some(0x01));
        assert_eq!(ht.get_base(&2), Some(0x02));
        assert_eq!(ht.get_base(&3), Some(0x03));
    }

    #[test]
    fn create_wide_geneaology() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);
        let _hts: Vec<HaplotypeRef> = (0..100)
            .map(|i| wt.create_descendant(vec![0], vec![Some(i)], 0))
            .collect();
        for (position, descendant) in wt.get_descendants().lock().iter().enumerate() {
            if let Some(d) = descendant.upgrade() {
                assert_eq!(d.get_base(&0), Some(position as u8));
            } else {
                panic!();
            }
        }
    }

    #[test]
    #[serial]
    fn single_recombination() {
        let mut haplotypes: Vec<HaplotypeRef> = Vec::new();
        let bytes = vec![Some(0x01); 100];
        let wildtype = Wildtype::new(bytes);
        haplotypes.push(wildtype.clone());
        for i in 0..100 {
            let ht = haplotypes.last().unwrap().clone();
            haplotypes.push(ht.create_descendant(vec![i], vec![Some(0x02)], 0));
        }
        haplotypes.push(wildtype.create_descendant(vec![0], vec![Some(0x03)], 0));
        for i in 1..100 {
            let ht = haplotypes.last().unwrap().clone();
            haplotypes.push(ht.create_descendant(vec![i], vec![Some(0x03)], 0));
        }
        let left_ancestor = haplotypes[100].clone();
        let right_ancestor = haplotypes[200].clone();
        let recombinant = Haplotype::create_recombinant(&left_ancestor, &right_ancestor, 10, 90, 0);

        let mut expected = vec![Some(0x03); 10];
        expected.append(&mut vec![Some(0x02); 80]);
        expected.append(&mut vec![Some(0x03); 10]);

        assert_eq!(recombinant.get_sequence(), expected);
    }

    #[test]
    #[serial]
    fn get_length() {
        let wildtype = Wildtype::new(vec![Some(0x00); 100]);
        let haplotype = wildtype.create_descendant(vec![0], vec![Some(0x01)], 0);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75, 0);
        assert_eq!(wildtype.get_length(), 100);
        assert_eq!(haplotype.get_length(), 100);
        assert_eq!(recombinant.get_length(), 100);
    }

    #[test]
    #[serial]
    fn get_sequence() {
        let wildtype = Wildtype::new(vec![Some(0x00); 100]);
        let haplotype = wildtype.create_descendant(vec![0], vec![Some(0x01)], 0);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75, 0);

        let mut expected = vec![Some(0x01)];
        expected.append(&mut vec![Some(0x00); 99]);

        assert_eq!(recombinant.get_sequence(), expected);
    }

    #[test]
    #[serial]
    fn create_tree() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);

        let ht = wt.create_descendant(vec![0], vec![Some(0x03)], 1);
        let ht_id = ht.get_block_id();
        assert_eq!(wt.get_tree(), format!("('{}':1)wt;", ht_id));

        let rc = Haplotype::create_recombinant(&wt, &ht, 1, 2, 2);
        let rc_id = rc.get_block_id();
        assert_eq!(
            wt.get_tree(),
            format!("((#R'{rc_id}':1)'{ht_id}':1,#R'{rc_id}':2)wt;")
        );
    }

    #[test]
    fn merge_nodes() {
        let bytes = vec![Some(0x00), Some(0x00), Some(0x00), Some(0x00)];
        let wt = Wildtype::new(bytes);

        let ht1 = wt.create_descendant(vec![0], vec![Some(0x01)], 1);
        let ht2 = ht1.create_descendant(vec![0], vec![Some(0x02)], 2);
        let ht3 = ht2.create_descendant(vec![0], vec![Some(0x03)], 3);

        let ht1weak = ht1.get_weak();
        let ht2weak = ht2.get_weak();

        // get ht1 from descendants as d1
        let d1 = wt
            .get_descendants()
            .lock()
            .first()
            .unwrap()
            .upgrade()
            .unwrap();

        assert_eq!(d1, ht1);

        // drop d1
        // this should not trigger nodes to be merged
        drop(d1);

        assert_eq!(ht2.try_get_ancestor().unwrap(), &ht1);

        // drop the outside references to ht1 and ht2
        // this should trigger ht1 to be merged with ht2 first
        // and then ht2 to be merged with ht3

        drop(ht1);

        assert!(!ht1weak.exists());

        drop(ht2);

        assert!(!ht2weak.exists());

        assert_eq!(ht3.get_generation(), 3);
        assert_eq!(ht3.get_mutations().len(), 1);
        assert_eq!(ht3.get_mutations().get(&0), Some(&Some(0x03)));
    }
}
