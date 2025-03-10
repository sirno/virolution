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
use parking_lot::{ReentrantMutex, ReentrantMutexGuard};
use seq_io::fasta::OwnedRecord;
use smallvec::SmallVec;
use std::cell::Cell;
use std::collections::HashMap;
use std::fmt;
use std::mem::ManuallyDrop;
use std::sync::atomic::{AtomicIsize, Ordering};
use std::sync::Arc;

use crate::encoding::Symbol;
use crate::errors::Result;
use crate::references::DescendantsCell;
use crate::references::{HaplotypeRef, HaplotypeRefTrait, HaplotypeWeak, HaplotypeWeakTrait};
use virolution_macros::require_deferred_drop;

use super::attributes::{AttributeSet, AttributeValue};
use super::cache::{CachedValue, VirolutionCache};
use super::AttributeSetDefinition;

// #[derive(Clone, Debug, Deref)]
// pub type Symbol = Option<u8>;

const SMALL_VEC_SIZE: usize = 1;
#[cfg(all(feature = "parallel", test))]
static N_LEAKED: AtomicIsize = AtomicIsize::new(0);
#[cfg(all(feature = "parallel", test))]
static N_DROPPED: AtomicIsize = AtomicIsize::new(0);

#[derive(Debug, Clone)]
pub struct Change<S: Symbol> {
    pub position: usize,
    pub from: S,
    pub to: S,
}

pub type Changes<S> = [Change<S>; SMALL_VEC_SIZE];

#[derive(Debug)]
pub enum Haplotype<S: Symbol> {
    Wildtype(Wildtype<S>),
    Mutant(Mutant<S>),
    Recombinant(Recombinant<S>),
}

#[derive(Debug)]
pub struct Wildtype<S: Symbol> {
    // head
    reference: HaplotypeWeak<S>,
    descendants: DescendantsCell<S>,

    // body
    sequence: Vec<S>,
    attributes: AttributeSet<S>,

    // sync
    // the number of descendants that have been dropped as a negative number
    _dirty_descendants: AtomicIsize,
}

/// Private struct to handle deferred drops
///
/// Contains a shareable lock to synchronize any deferred drops and a cell to temporarily leak
/// self-reference to the haplotype
#[derive(Derivative)]
#[derivative(Debug)]
struct DeferenceCell<T> {
    // a lock to synchronize deferred drops, needs to be shared between the original and the
    // deferred haplotype to ensure that the deference requesters can both, access through the
    // reference and the deferred haplotype directly
    lock: Arc<ReentrantMutex<Cell<u8>>>,
    dirty: Cell<bool>,
    #[derivative(Debug = "ignore")]
    leak: Cell<Option<ManuallyDrop<T>>>,
}

impl<T> DeferenceCell<T> {
    fn new() -> Self {
        Self {
            lock: Arc::new(ReentrantMutex::new(Cell::new(0))),
            leak: Cell::new(None),
            dirty: Cell::new(false),
        }
    }

    fn share(&self) -> Self {
        Self {
            lock: self.lock.clone(),
            leak: Cell::new(None),
            dirty: Cell::new(false),
        }
    }

    fn lock(&self) -> ReentrantMutexGuard<Cell<u8>> {
        self.lock.lock()
    }

    fn require_deferred_drop(&self) {
        let guard = self.lock();
        guard.set(guard.get() + 1);
    }

    // returns true if a deferred drop should be executed
    fn inquire_deferred_drop(&self) -> bool {
        let guard = self.lock();

        // decrement the number of requests
        guard.set(guard.get() - 1);

        // return true if the haplotype is dirty and no more requests are pending
        self.dirty.get() && guard.get() == 0
    }

    fn get_requests(&self) -> u8 {
        self.lock.lock().get()
    }

    fn is_dirty(&self) -> bool {
        self.dirty.get()
    }

    fn set_dirty(&self) {
        self.dirty.set(true);
    }
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Mutant<S: Symbol> {
    // head
    reference: HaplotypeWeak<S>,
    wildtype: HaplotypeWeak<S>,
    ancestor: HaplotypeRef<S>,
    descendants: DescendantsCell<S>,

    // body
    changes: SmallVec<Changes<S>>,
    attributes: AttributeSet<S>,

    // sync
    // the number of descendants that have been dropped as a negative number
    _dirty_descendants: AtomicIsize,
    _defer_drop: DeferenceCell<HaplotypeRef<S>>,
    // synchronization for merging nodes while allowing for parallel access
    // this field will be used to request deferred drops, when it is non-zero, merges will not
    // drop but instead defer the drop to the setter
    // any thread that requires a deferred drop is responsible for decrementing this field before
    // dropping the reference (safe handling is ensured by the `require_deferred_drop` attribute)
    // _defer_drop: Arc<Mutex<u8>>,
    // this field will be used to create a self-reference to the haplotype after it has been
    // removed from the tree, this allows us to safely drop the haplotype when it is no longer used
    // by consuming the reference when it is no longer used
    // #[derivative(Debug = "ignore")]
    // _drop: Cell<Option<HaplotypeRef<S>>>,
}

#[derive(Debug)]
pub struct Recombinant<S: Symbol> {
    // head
    reference: HaplotypeWeak<S>,
    wildtype: HaplotypeWeak<S>,
    left_ancestor: HaplotypeRef<S>,
    right_ancestor: HaplotypeRef<S>,
    descendants: DescendantsCell<S>,

    // body
    left_position: usize,
    right_position: usize,
    attributes: AttributeSet<S>,

    // sync
    // the number of descendants that have been dropped as a negative number
    _dirty_descendants: AtomicIsize,
}

// We have to explicitly implement Send and Sync for Haplotype, because it contains several
// read-only fields that are not thread-safe. For example, it may contain a HashMap, which reflects
// the changes in the haplotype data, that is created during the construction of the Haplotype.
// To my knowledge the only effect of these lines (right now) is that it silences the warning
// that Arc references an object that is not Send and Sync.
#[cfg(feature = "parallel")]
unsafe impl<S: Symbol> Send for Haplotype<S> {}
#[cfg(feature = "parallel")]
unsafe impl<S: Symbol> Sync for Haplotype<S> {}

impl<S: Symbol> Drop for Haplotype<S> {
    fn drop(&mut self) {
        match self {
            Haplotype::Wildtype(_wt) => {}
            Haplotype::Mutant(mt) => {
                if !mt._defer_drop.is_dirty() {
                    mt.ancestor.add_dirty_descendants();
                }
            }
            Haplotype::Recombinant(rc) => {
                rc.left_ancestor.add_dirty_descendants();
                rc.right_ancestor.add_dirty_descendants();
            }
        }
    }
}

impl<S: Symbol> Haplotype<S> {
    /// Create a new haplotype from a sequence and an attribute definition.
    ///
    /// This can be the starting point of a new phylogeny. The attribute definition is used to
    /// register different traits that can be used to derive new attributes for the haplotype.
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        sequence: Vec<S>,
        attribute_definition: &AttributeSetDefinition<S>,
    ) -> HaplotypeRef<S> {
        Wildtype::new(sequence, attribute_definition)
    }

    /// Create a mutant haplotype from a parent haplotype and a set of changes.
    ///
    /// The changes are represented as a vector of positions and a vector of symbols. The parent
    /// haplotype is the ancestor of the new mutant haplotype.
    pub fn create_descendant(&self, positions: Vec<usize>, changes: Vec<S>) -> HaplotypeRef<S> {
        let ancestor = self.get_reference();
        let wildtype = self.get_wildtype();

        let changes = positions
            .iter()
            .zip(changes.iter())
            .map(|(&position, &to)| {
                let from = self.get_base(&position);
                assert_ne!(from, to);
                Change { position, from, to }
            })
            .collect();

        let descendant = Mutant::new(ancestor, wildtype, changes);

        self.add_descendant(descendant.get_weak());

        descendant
    }

    /// Create a recombinant haplotype from two parent haplotypes and two positions.
    ///
    /// The two parent haplotypes are the ancestors of the new recombinant haplotype. The
    /// recombinant is characterized by two positions, that characterize the left and right
    /// starting points of the recombination, i.e. the left parent haplotype will be used from the
    /// left position to the right position and the right parent haplotype will be used from the
    /// right position to the left position.
    pub fn create_recombinant(
        left_ancestor: &HaplotypeRef<S>,
        right_ancestor: &HaplotypeRef<S>,
        left_position: usize,
        right_position: usize,
    ) -> HaplotypeRef<S> {
        let wildtype = left_ancestor.get_wildtype();

        let recombinant = Recombinant::new(
            wildtype,
            left_ancestor.clone(),
            right_ancestor.clone(),
            left_position,
            right_position,
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
    pub(crate) fn try_unwrap_mutant(&self) -> Option<&Mutant<S>> {
        match self {
            Haplotype::Mutant(ht) => Some(ht),
            _ => None,
        }
    }

    /// Returns a reference to the changes that are present in the haplotype if the type allows it.
    pub fn try_get_changes(&self) -> Option<&SmallVec<Changes<S>>> {
        match self {
            Haplotype::Mutant(ht) => Some(&ht.changes),
            _ => None,
        }
    }

    /// Returns a reference to an ancestor if the type allows it.
    ///
    /// If there are multiple ancestors, in the future any ancestor may be returned.
    pub fn try_get_ancestor(&self) -> Option<&HaplotypeRef<S>> {
        match self {
            Haplotype::Mutant(ht) => Some(&ht.ancestor),
            Haplotype::Recombinant(rc) => Some(&rc.left_ancestor),
            _ => None,
        }
    }

    pub fn get_reference(&self) -> HaplotypeRef<S> {
        let weak = match self {
            Haplotype::Wildtype(wt) => &wt.reference,
            Haplotype::Mutant(ht) => &ht.reference,
            Haplotype::Recombinant(rc) => &rc.reference,
        };
        weak.upgrade().expect("Self-reference has been dropped.")
    }

    pub fn get_wildtype(&self) -> HaplotypeWeak<S> {
        match self {
            Haplotype::Wildtype(wt) => wt.get_reference().get_weak(),
            Haplotype::Mutant(ht) => ht.wildtype.clone(),
            Haplotype::Recombinant(rc) => rc.wildtype.clone(),
        }
    }

    pub fn get_ancestors(&self) -> (Option<HaplotypeRef<S>>, Option<HaplotypeRef<S>>) {
        match self {
            Haplotype::Mutant(ht) => (Some(ht.ancestor.clone()), None),
            Haplotype::Recombinant(rc) => (
                Some(rc.left_ancestor.clone()),
                Some(rc.right_ancestor.clone()),
            ),
            _ => (None, None),
        }
    }

    pub fn get_descendants(&self) -> &DescendantsCell<S> {
        match self {
            Haplotype::Wildtype(wt) => &wt.descendants,
            Haplotype::Mutant(ht) => &ht.descendants,
            Haplotype::Recombinant(rc) => &rc.descendants,
        }
    }

    /// Returns the number of descendants.
    fn add_descendant(&self, descendant: HaplotypeWeak<S>) {
        let (descendants, dirty_descendants) = match self {
            Haplotype::Wildtype(wt) => (&wt.descendants, &wt._dirty_descendants),
            Haplotype::Mutant(ht) => (&ht.descendants, &ht._dirty_descendants),
            Haplotype::Recombinant(rc) => (&rc.descendants, &rc._dirty_descendants),
        };

        let mut descendants_guard = descendants.lock();

        // if there are dirty descendants, replace one of them
        if dirty_descendants.load(Ordering::Acquire) < 0
            && let Some(idx) = descendants_guard.iter().position(|x| !x.exists())
        {
            descendants_guard[idx] = descendant;
            dirty_descendants.fetch_add(1, Ordering::AcqRel);
            return;
        }

        descendants_guard.push(descendant);
    }

    fn add_dirty_descendants(&self) {
        let dirty_descendants = match self {
            Haplotype::Wildtype(wt) => &wt._dirty_descendants,
            Haplotype::Mutant(ht) => &ht._dirty_descendants,
            Haplotype::Recombinant(rc) => &rc._dirty_descendants,
        };
        dirty_descendants.fetch_sub(1, Ordering::AcqRel);
    }

    pub fn get_base(&self, position: &usize) -> S {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Mutant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_sequence(&self) -> Vec<S> {
        let mutations = self.get_mutations();
        let mut sequence = self.get_wildtype_sequence();

        for (&position, &new) in mutations.iter() {
            sequence[position] = S::decode(&new);
        }

        sequence
    }

    pub fn get_wildtype_sequence(&self) -> Vec<S> {
        match self {
            Haplotype::Wildtype(wt) => wt.sequence.to_vec(),
            Haplotype::Mutant(_ht) => self
                .get_wildtype()
                .upgrade()
                .unwrap()
                .get_wildtype_sequence(),
            Haplotype::Recombinant(_rc) => self
                .get_wildtype()
                .upgrade()
                .unwrap()
                .get_wildtype_sequence(),
        }
    }

    pub fn get_string(&self) -> String {
        let mutations = self.get_mutations();
        let wildtype = self.get_wildtype_sequence();

        let mut out = String::new();
        for (position, to) in mutations.iter() {
            let from = wildtype[*position];
            out.push_str(format!(";{position}:{from}->{}", S::decode(to)).as_str());
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
    pub fn get_mutations(&self) -> CachedValue<HashMap<usize, u8>> {
        match self {
            Haplotype::Wildtype(_wt) => Default::default(),
            Haplotype::Mutant(ht) => ht.get_mutations(),
            Haplotype::Recombinant(rc) => rc.get_mutations(),
        }
    }

    pub fn get_attributes(&self) -> &AttributeSet<S> {
        match self {
            Haplotype::Wildtype(wt) => &wt.attributes,
            Haplotype::Mutant(ht) => &ht.attributes,
            Haplotype::Recombinant(rc) => &rc.attributes,
        }
    }

    pub fn get_attribute(&self, id: &str) -> Option<AttributeValue> {
        self.get_attributes().get(id)
    }

    pub fn get_or_compute_attribute(&self, id: &'static str) -> Result<AttributeValue> {
        self.get_attributes().get_or_compute(id)
    }

    pub fn get_record(&self, head: &str) -> OwnedRecord {
        OwnedRecord {
            head: head.to_string().as_bytes().to_vec(),
            seq: self
                .get_sequence()
                .into_iter()
                .map(|symbol| symbol.encode())
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

    pub fn get_subtree(&self, ancestor: HaplotypeWeak<S>) -> String {
        match self {
            Haplotype::Wildtype(wt) => wt.get_subtree(),
            Haplotype::Mutant(ht) => ht.get_subtree(),
            Haplotype::Recombinant(rc) => rc.get_subtree(ancestor),
        }
    }
}

impl<S: Symbol> Wildtype<S> {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        sequence: Vec<S>,
        attribute_definition: &AttributeSetDefinition<S>,
    ) -> HaplotypeRef<S> {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Wildtype(Self {
                // head
                reference: reference.clone(),
                descendants: DescendantsCell::new(),

                // body
                sequence: sequence.clone(),
                attributes: attribute_definition.create(reference.clone()),

                // sync
                _dirty_descendants: AtomicIsize::new(0),
            })
        })
    }

    pub fn get_reference(&self) -> HaplotypeRef<S> {
        self.reference
            .upgrade()
            .expect("Self-reference has been dropped.")
    }

    pub fn get_base(&self, position: &usize) -> S {
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

impl<S: Symbol> Mutant<S> {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        ancestor: HaplotypeRef<S>,
        wildtype: HaplotypeWeak<S>,
        changes: SmallVec<Changes<S>>,
    ) -> HaplotypeRef<S> {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Mutant(Self {
                // head
                reference: reference.clone(),
                wildtype: wildtype.clone(),
                ancestor: ancestor.clone(),
                descendants: DescendantsCell::new(),

                // body
                changes,
                attributes: ancestor.get_attributes().derive(reference.clone()),

                // sync
                _dirty_descendants: AtomicIsize::new(0),
                _defer_drop: DeferenceCell::new(),
            })
        })
    }

    /// Creates a new mutant haplotype and replaces `last` with it.
    ///
    /// This function is not thread-safe and should only be called by a single thread at a time.
    pub(crate) unsafe fn new_and_replace(
        last: &Self,
        ancestor: HaplotypeRef<S>,
        wildtype: HaplotypeWeak<S>,
        changes: SmallVec<Changes<S>>,
    ) {
        // collect all descendants that are still alive
        let descendants: Vec<HaplotypeWeak<S>> = last
            .descendants
            .lock()
            .iter()
            .filter(|x| x.exists())
            .cloned()
            .collect();

        // make the ancestor aware of the new descendant
        ancestor.add_descendant(last.reference.clone());

        // create new node
        let tmp_ref = HaplotypeRef::new(Haplotype::Mutant(Self {
            // head
            reference: last.reference.clone(),
            wildtype,
            ancestor,
            descendants: DescendantsCell::from_iter(descendants),

            // body
            changes,
            attributes: last.attributes.clone(),

            // sync
            _dirty_descendants: AtomicIsize::new(0),
            _defer_drop: last._defer_drop.share(),
            // _drop: Cell::new(None),
        }));

        let old_ptr = last.reference.as_ptr() as *mut Mutant<S>;
        let new_ptr = tmp_ref.as_ptr() as *mut Mutant<S>;

        // synchronize swapping and deference with other threads
        let _guard = (*old_ptr)._defer_drop.lock();
        let _other_guard = (*new_ptr)._defer_drop.lock();

        // replace the old reference with the new one inside the haplotype ref
        // this swaps the reference count of the old reference for the new one
        std::ptr::swap(old_ptr, new_ptr);

        // mark the old reference as dirty
        (*old_ptr)._defer_drop.set_dirty();

        // check if anyone has requested a deferred drop
        if (*old_ptr)._defer_drop.get_requests() > 0 {
            eprintln!("Leaking reference to haplotype.");
            let dangled_ref = ManuallyDrop::new(tmp_ref);
            (*old_ptr)._defer_drop.leak.set(Some(dangled_ref));
            #[cfg(all(feature = "parallel", test))]
            N_LEAKED.fetch_add(1, Ordering::Relaxed);
        }
    }

    #[require_deferred_drop]
    pub fn get_base(&self, position: &usize) -> S {
        match self.changes.iter().find(|x| x.position == *position) {
            Some(change) => change.to,
            None => self.ancestor.get_base(position),
        }
    }

    pub fn iter_changes(&self) -> impl Iterator<Item = &Change<S>> + '_ {
        self.changes.iter()
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    pub fn get_ancestor(&self) -> HaplotypeRef<S> {
        self.ancestor.clone()
    }

    pub fn get_changes(&self) -> &SmallVec<Changes<S>> {
        &self.changes
    }

    /// Returns a HashMap of mutations that are present in the haplotype.
    #[require_deferred_drop]
    pub fn get_mutations(&self) -> CachedValue<HashMap<usize, u8>> {
        let key = &self.changes as *const _;
        let mut mutations = self.ancestor.get_mutations().clone_inner();

        let wt_ref = self.wildtype.upgrade().unwrap();

        assert_eq!(&self.changes as *const _, key);

        self.changes.iter().for_each(|change| {
            let wt_base = wt_ref.get_base(&change.position);
            if change.to == wt_base {
                mutations.remove(&change.position);
            } else {
                mutations.insert(change.position, *change.to.index() as u8);
            }
        });

        CachedValue::new(mutations)
    }

    pub fn get_subtree(&self) -> String {
        let current = self.reference.upgrade().unwrap();

        let block_id = current.get_block_id();

        let node = match (
            self.attributes.get("generation"),
            self.ancestor.get_attributes().get("generation"),
        ) {
            (
                Some(AttributeValue::U64(generation)),
                Some(AttributeValue::U64(ancestor_generation)),
            ) => {
                let branch_length = generation - ancestor_generation;
                format!("'{block_id}':{branch_length}")
            }
            _ => format!("'{block_id}'"),
        };

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
            .checked_add_signed(self._dirty_descendants.load(Ordering::Acquire))
            .expect("Number of descendants overflowed.");

        let descendant = match n_descendants {
            1 => self
                .descendants
                .lock()
                .iter()
                .filter_map(|x| x.upgrade())
                .next(),
            _ => return,
        };

        let descendant_inner = match descendant.as_ref().and_then(|x| x.try_unwrap_mutant()) {
            Some(x) => x,
            None => return,
        };

        let merger: [&Mutant<S>; 2] = [descendant_inner, self];

        // aggregate changes
        let changes: SmallVec<Changes<S>> = merger
            .iter()
            .rev()
            .flat_map(|x| x.changes.iter())
            .cloned()
            .collect();

        // extract ancestor and wildtype
        let ancestor = self.ancestor.clone();
        let wildtype = self.wildtype.clone();

        // create new node
        unsafe {
            Mutant::new_and_replace(descendant_inner, ancestor, wildtype, changes);
        };
    }

    /// Notifies that any drop needs to be deferred
    fn require_deferred_drop(&self) {
        self._defer_drop.require_deferred_drop();
    }

    /// Check if a drop has been deferred and if any other thread has requested deferred drop.
    /// If no other thread has requested deferred drop, the drop will be executed.
    fn inquire_deferred_drop(&self) {
        if self._defer_drop.inquire_deferred_drop() {
            if let Some(mut dangling_ptr) = self._defer_drop.leak.take() {
                #[cfg(all(feature = "parallel", test))]
                N_DROPPED.fetch_add(1, Ordering::Relaxed);
                eprintln!("Dropping reference to haplotype.");
                unsafe {
                    ManuallyDrop::drop(&mut dangling_ptr);
                }
            }
        }
    }
}

impl<S: Symbol> Recombinant<S> {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        wildtype: HaplotypeWeak<S>,
        left_ancestor: HaplotypeRef<S>,
        right_ancestor: HaplotypeRef<S>,
        left_position: usize,
        right_position: usize,
    ) -> HaplotypeRef<S> {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Recombinant(Self {
                // head
                reference: reference.clone(),
                wildtype: wildtype.clone(),
                left_ancestor: left_ancestor.clone(),
                right_ancestor: right_ancestor.clone(),
                descendants: DescendantsCell::new(),

                // body
                left_position,
                right_position,
                attributes: left_ancestor.get_attributes().derive(reference.clone()),

                // sync
                _dirty_descendants: AtomicIsize::new(0),
            })
        })
    }

    pub fn get_base(&self, position: &usize) -> S {
        if *position >= self.left_position && *position < self.right_position {
            return self.left_ancestor.get_base(position);
        }

        self.right_ancestor.get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    fn get_mutations_cache(
    ) -> &'static ::cached::once_cell::sync::Lazy<VirolutionCache<HashMap<usize, u8>>> {
        // at this point, we cannot use generics in statics, so we encode the data for now until
        // later or const generics become available
        static MUTATIONS_CACHE: ::cached::once_cell::sync::Lazy<
            VirolutionCache<HashMap<usize, u8>>,
        > = ::cached::once_cell::sync::Lazy::new(|| VirolutionCache::new(100));
        &MUTATIONS_CACHE
    }

    pub fn get_mutations(&self) -> CachedValue<HashMap<usize, u8>> {
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

    pub fn get_subtree(&self, ancestor: HaplotypeWeak<S>) -> String {
        let block_id = self.reference.get_block_id();

        let node = match (
            self.attributes.get("generation"),
            ancestor
                .upgrade()
                .unwrap()
                .get_attributes()
                .get("generation"),
        ) {
            (
                Some(AttributeValue::U64(generation)),
                Some(AttributeValue::U64(ancestor_generation)),
            ) => {
                let branch_length = generation - ancestor_generation;
                format!("#R'{block_id}':{branch_length}")
            }
            _ => format!("#R'{block_id}'"),
        };

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

impl<S: Symbol> fmt::Display for Haplotype<S> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.get_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::attributes::AttributeSetDefinition;
    use crate::encoding::Nucleotide as Nt;
    use crate::providers::Generation;
    use serial_test::serial;
    use std::sync::Arc;
    #[cfg(feature = "parallel")]
    use std::sync::Mutex;

    #[test]
    fn initiate_wildtype() {
        let symbols = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(symbols.clone(), &attribute_definition);
        for i in symbols.iter().enumerate() {
            assert_eq!(wt.get_base(&i.0), *i.1);
        }
    }

    #[test]
    fn create_descendant() {
        let symbols = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(symbols.clone(), &attribute_definition);
        let ht = wt.create_descendant(vec![0], vec![Nt::G]);
        assert_eq!(ht.get_base(&0), Nt::G);
        assert_eq!(ht.get_base(&1), Nt::T);
        assert_eq!(ht.get_base(&2), Nt::C);
        assert_eq!(ht.get_base(&3), Nt::G);
    }

    #[test]
    fn create_wide_geneaology() {
        let symbols = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(symbols.clone(), &attribute_definition);
        let _hts: Vec<HaplotypeRef<Nt>> = (0..100)
            .map(|i| {
                wt.create_descendant(
                    vec![0],
                    vec![Nt::decode(&(((i % (Nt::SIZE - 1)) + 1) as u8))],
                )
            })
            .collect();
        for (i, descendant) in wt.get_descendants().lock().iter().enumerate() {
            if let Some(d) = descendant.upgrade() {
                assert_eq!(
                    d.get_base(&0),
                    Nt::decode(&(((i % (Nt::SIZE - 1)) + 1) as u8))
                );
            } else {
                panic!();
            }
        }
    }

    #[test]
    #[serial]
    fn single_recombination() {
        let mut haplotypes: Vec<HaplotypeRef<Nt>> = Vec::new();
        let symbols = vec![Nt::T; 100];
        let attribute_definition = AttributeSetDefinition::new();
        let wildtype = Wildtype::new(symbols, &attribute_definition);
        haplotypes.push(wildtype.clone());
        for i in 0..100 {
            let ht = haplotypes.last().unwrap().clone();
            haplotypes.push(ht.create_descendant(vec![i], vec![Nt::C]));
        }
        haplotypes.push(wildtype.create_descendant(vec![0], vec![Nt::G]));
        for i in 1..100 {
            let ht = haplotypes.last().unwrap().clone();
            haplotypes.push(ht.create_descendant(vec![i], vec![Nt::G]));
        }
        let left_ancestor = haplotypes[100].clone();
        let right_ancestor = haplotypes[200].clone();
        let recombinant = Haplotype::create_recombinant(&left_ancestor, &right_ancestor, 10, 90);

        let mut expected = vec![Nt::G; 10];
        expected.append(&mut vec![Nt::C; 80]);
        expected.append(&mut vec![Nt::G; 10]);

        assert_eq!(recombinant.get_sequence(), expected);
    }

    #[test]
    #[serial]
    fn get_length() {
        let symbols = vec![Nt::A; 100];
        let attribute_definition = AttributeSetDefinition::new();
        let wildtype = Wildtype::new(symbols, &attribute_definition);
        let haplotype = wildtype.create_descendant(vec![0], vec![Nt::T]);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75);
        assert_eq!(wildtype.get_length(), 100);
        assert_eq!(haplotype.get_length(), 100);
        assert_eq!(recombinant.get_length(), 100);
    }

    #[test]
    #[serial]
    fn get_sequence() {
        let symbols = vec![Nt::A; 100];
        let attribute_definition = AttributeSetDefinition::new();
        let wildtype = Wildtype::new(symbols, &attribute_definition);
        let haplotype = wildtype.create_descendant(vec![0], vec![Nt::T]);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75);

        let mut expected = vec![Nt::T];
        expected.append(&mut vec![Nt::A; 99]);

        assert_eq!(recombinant.get_sequence(), expected);
    }

    #[test]
    #[serial]
    fn create_tree() {
        let symbols = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let mut attribute_definition = AttributeSetDefinition::new();
        let generation_provider = Arc::new(Generation::new(0));
        attribute_definition.register(generation_provider.clone());
        let wt = Wildtype::new(symbols, &attribute_definition);

        generation_provider.increment();
        let ht = wt.create_descendant(vec![0], vec![Nt::G]);
        let ht_id = ht.get_block_id();
        assert_eq!(wt.get_tree(), format!("('{}':1)wt;", ht_id));

        generation_provider.increment();
        let rc = Haplotype::create_recombinant(&wt, &ht, 1, 2);
        let rc_id = rc.get_block_id();
        assert_eq!(
            wt.get_tree(),
            format!("((#R'{rc_id}':1)'{ht_id}':1,#R'{rc_id}':2)wt;")
        );
    }

    #[test]
    fn merge_nodes() {
        let symbols = vec![Nt::A, Nt::T, Nt::C, Nt::G];
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(symbols, &attribute_definition);

        let ht1 = wt.create_descendant(vec![0], vec![Nt::T]);
        let ht2 = ht1.create_descendant(vec![0], vec![Nt::C]);
        let ht3 = ht2.create_descendant(vec![0], vec![Nt::G]);

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

        assert_eq!(ht3.get_mutations().len(), 1);
        assert_eq!(ht3.get_mutations().get(&0), Some(&(*Nt::G.index() as u8)));
    }

    #[cfg(feature = "parallel")]
    fn mutator(
        pop: &[Mutex<HaplotypeRef<Nt>>],
        pop_size: usize,
        n_mutations: usize,
        n_sites: usize,
        n_symbols: usize,
    ) {
        use rand::prelude::*;

        let mut rng = rand::rng();

        for i in 0..n_mutations {
            let from = rng.random_range(0..pop_size);
            let to = rng.random_range(0..pop_size);

            let descendant = {
                let ht = pop[from].lock().expect("Failed to lock ancestor.");
                let pos = i % n_sites;
                let sym = ht.get_base(&pos);
                ht.create_descendant(
                    vec![pos],
                    vec![Nt::decode(&(((sym.index() + 1) % n_symbols) as u8))],
                )
            };

            let mut field = pop[to].lock().expect("Failed to lock field to change.");
            *field = descendant;
        }
    }

    #[cfg(feature = "parallel")]
    fn reader(pop: &[Mutex<HaplotypeRef<Nt>>], pop_size: usize, n_reads: usize, n_sites: usize) {
        use rand::prelude::*;
        for _ in 0..n_reads {
            let mut rng = rand::rng();
            let index = rng.random_range(0..pop_size);
            let ht = pop[index].lock().expect("Failed to lock haplotype.");
            let sequence = ht.get_sequence();
            assert_eq!(sequence.len(), n_sites);
        }
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn merge_nodes_stress() {
        use std::sync::Mutex;
        use std::thread;

        let n_sites = 7;
        let n_symbols = 4;

        let symbols = vec![Nt::A; n_sites];
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(symbols, &attribute_definition);

        let n_mutations = 100000;
        let n_reads = 10000;

        let pop_size = 11;
        let pop: Vec<Mutex<HaplotypeRef<Nt>>> =
            (0..pop_size).map(|_| Mutex::new(wt.clone())).collect();

        thread::scope(|s| {
            s.spawn(|| mutator(&pop, pop_size, n_mutations, n_sites, n_symbols));
            s.spawn(|| reader(&pop, pop_size, n_reads, n_sites));
        });

        thread::sleep(std::time::Duration::from_millis(100));

        let n_leaked = N_LEAKED.load(Ordering::Acquire);
        let n_dropped = N_DROPPED.load(Ordering::Acquire);
        assert_eq!(
            n_leaked, n_dropped,
            "Number of leaked and dropped nodes should be equal."
        );
    }

    #[test]
    #[cfg(feature = "parallel")]
    fn merge_nodes_high_stress() {
        use std::sync::Mutex;
        use std::thread;

        // TODO: investigate potential smallvec issue (unreachable code entered?)

        let n_sites = 7;
        let n_symbols = 4;

        let symbols = vec![Nt::A; n_sites];
        let attribute_definition = AttributeSetDefinition::new();
        let wt = Wildtype::new(symbols, &attribute_definition);

        let n_mutations = 1000;
        let n_reads = 100;

        let pop_size = 100;
        let pop: Vec<Mutex<HaplotypeRef<Nt>>> =
            (0..pop_size).map(|_| Mutex::new(wt.clone())).collect();

        thread::scope(|s| {
            s.spawn(|| mutator(&pop, pop_size, n_mutations, n_sites, n_symbols));
            // s.spawn(|| mutator(&pop, pop_size, n_mutations, n_sites, n_symbols));
            // s.spawn(|| mutator(&pop, pop_size, n_mutations, n_sites, n_symbols));
            // s.spawn(|| mutator(&pop, pop_size, n_mutations, n_sites, n_symbols));
            s.spawn(|| reader(&pop, pop_size, n_reads, n_sites));
            // s.spawn(|| reader(&pop, pop_size, n_reads, n_sites));
            // s.spawn(|| reader(&pop, pop_size, n_reads, n_sites));
            // s.spawn(|| reader(&pop, pop_size, n_reads, n_sites));
        });
    }
}
