use super::fitness::FitnessTable;
use super::references::DescendantsCell;
use super::references::{HaplotypeRef, HaplotypeWeak};
use derivative::Derivative;
use phf::phf_map;
use seq_io::fasta::OwnedRecord;
use std::collections::HashMap;
use std::fmt;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::OnceLock;

// #[derive(Clone, Debug, Deref)]
pub type Symbol = Option<u8>;

pub static FASTA_ENCODE: phf::Map<u8, u8> = phf_map! {
    0x00u8 => 0x41,
    0x01u8 => 0x47,
    0x02u8 => 0x54,
    0x03u8 => 0x43,
};
pub static FASTA_DECODE: phf::Map<u8, u8> = phf_map! {
    0x41u8 => 0x00,
    0x47u8 => 0x01,
    0x54u8 => 0x02,
    0x43u8 => 0x03,
};

#[derive(Debug)]
pub enum Haplotype {
    Wildtype(Wildtype),
    Descendant(Descendant),
    Recombinant(Recombinant),
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Wildtype {
    reference: HaplotypeWeak,
    sequence: Vec<Symbol>,
    descendants: DescendantsCell,
    dirty_descendants: AtomicUsize,
}

#[derive(Debug)]
pub struct Descendant {
    reference: HaplotypeWeak,
    wildtype: HaplotypeWeak,
    ancestor: HaplotypeRef,
    changes: HashMap<usize, (Symbol, Symbol)>,
    generation: usize,
    fitness: OnceLock<f64>,
    descendants: DescendantsCell,
    dirty_descendants: AtomicUsize,
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
    fitness: OnceLock<f64>,
    descendants: DescendantsCell,
    dirty_descendants: AtomicUsize,
}

impl Drop for Haplotype {
    fn drop(&mut self) {
        match self {
            Haplotype::Wildtype(_wt) => {}
            Haplotype::Descendant(ht) => ht.ancestor.increment_dirty_descendants(),
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

        let descendant = Descendant::new(ancestor, wildtype.get_weak(), changes, generation);

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

    pub fn is_wildtype(&self) -> bool {
        matches!(self, Haplotype::Wildtype(_))
    }

    pub fn is_descendant(&self) -> bool {
        matches!(self, Haplotype::Descendant(_))
    }

    pub fn is_recombinant(&self) -> bool {
        matches!(self, Haplotype::Recombinant(_))
    }

    fn unwrap_mutant(&self) -> &Descendant {
        match self {
            Haplotype::Descendant(ht) => ht,
            _ => panic!("Haplotype is not a mutant."),
        }
    }

    fn try_unwrap_mutant(&self) -> Option<&Descendant> {
        match self {
            Haplotype::Descendant(ht) => Some(ht),
            _ => None,
        }
    }

    pub fn get_reference(&self) -> HaplotypeRef {
        let weak = match self {
            Haplotype::Wildtype(wt) => &wt.reference,
            Haplotype::Descendant(ht) => &ht.reference,
            Haplotype::Recombinant(rc) => &rc.reference,
        };
        weak.upgrade().expect("Self-reference has been dropped.")
    }

    pub fn get_wildtype(&self) -> HaplotypeRef {
        match self {
            Haplotype::Wildtype(wt) => wt.get_reference().clone(),
            Haplotype::Descendant(ht) => ht.wildtype.upgrade().unwrap().clone(),
            Haplotype::Recombinant(rc) => rc.wildtype.upgrade().unwrap().clone(),
        }
    }

    fn get_ancestors(&self) -> (Option<HaplotypeRef>, Option<HaplotypeRef>) {
        match self {
            Haplotype::Descendant(ht) => (Some(ht.ancestor.clone()), None),
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
            Haplotype::Descendant(ht) => &ht.descendants,
            Haplotype::Recombinant(rc) => &rc.descendants,
        }
    }

    pub fn n_descendants(&self) -> usize {
        match self {
            Haplotype::Wildtype(wt) => {
                wt.descendants.lock().len() - wt.dirty_descendants.load(Ordering::Relaxed)
            }
            Haplotype::Descendant(ht) => {
                ht.descendants.lock().len() - ht.dirty_descendants.load(Ordering::Relaxed)
            }
            Haplotype::Recombinant(rc) => {
                rc.descendants.lock().len() - rc.dirty_descendants.load(Ordering::Relaxed)
            }
        }
    }

    fn add_descendant(&self, descendant: HaplotypeWeak) {
        let (descendants, dirty_descendants) = match self {
            Haplotype::Wildtype(wt) => (&wt.descendants, &wt.dirty_descendants),
            Haplotype::Descendant(ht) => (&ht.descendants, &ht.dirty_descendants),
            Haplotype::Recombinant(rc) => (&rc.descendants, &rc.dirty_descendants),
        };

        let mut descendants_guard = descendants.lock();

        // if there are dirty descendants, replace one of them
        if dirty_descendants.load(Ordering::Relaxed) > 0 && let Some(idx) = descendants_guard
            .iter()
            .position(|x| !x.exists())
        {
            descendants_guard[idx] = descendant;
            dirty_descendants.fetch_sub(1, Ordering::Relaxed);
            return;
        }

        descendants_guard.push(descendant);
    }

    fn increment_dirty_descendants(&self) {
        let dirty_descendants = match self {
            Haplotype::Wildtype(wt) => &wt.dirty_descendants,
            Haplotype::Descendant(ht) => &ht.dirty_descendants,
            Haplotype::Recombinant(rc) => &rc.dirty_descendants,
        };
        dirty_descendants.fetch_add(1, Ordering::Relaxed);
    }

    pub fn get_base(&self, position: &usize) -> Symbol {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Descendant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_generation(&self) -> usize {
        match self {
            Haplotype::Wildtype(_wt) => 0,
            Haplotype::Descendant(ht) => ht.get_generation(),
            Haplotype::Recombinant(rc) => rc.get_generation(),
        }
    }

    pub fn get_sequence(&self) -> Vec<Symbol> {
        let mutations = self.get_mutations();
        let mut sequence = self.get_wildtype_sequence();

        for (position, (_, to)) in mutations {
            sequence[position] = to;
        }

        sequence
    }

    pub fn get_wildtype_sequence(&self) -> Vec<Symbol> {
        match self {
            Haplotype::Wildtype(wt) => wt.sequence.to_vec(),
            Haplotype::Descendant(_ht) => self.get_wildtype().get_wildtype_sequence(),
            Haplotype::Recombinant(_rc) => self.get_wildtype().get_wildtype_sequence(),
        }
    }

    pub fn get_string(&self) -> String {
        let mutations = self.get_mutations();
        let wildtype = self.get_wildtype_sequence();

        let mut out = String::new();
        for (position, (from, to)) in mutations.iter() {
            match (from, to) {
                (Some(f), Some(t)) => {
                    out.push_str(format!(";{position}:{f}->{t}").as_str());
                }
                (None, Some(t)) => {
                    if let Some(wt_symbol) = wildtype[*position] {
                        out.push_str(format!(";{position}:{wt_symbol}->{t}").as_str())
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

    pub fn get_mutations(&self) -> HashMap<usize, (Symbol, Symbol)> {
        match self {
            Haplotype::Wildtype(_wt) => HashMap::new(),
            Haplotype::Descendant(ht) => ht.get_mutations(),
            Haplotype::Recombinant(rc) => rc.get_mutations(),
        }
    }

    pub fn get_fitness(&self, fitness_table: &FitnessTable) -> f64 {
        match self {
            Haplotype::Wildtype(_wt) => 1.,
            Haplotype::Descendant(ht) => *ht.get_fitness(fitness_table),
            Haplotype::Recombinant(rc) => *rc.get_fitness(fitness_table),
        }
    }

    pub fn get_record(&self, head: &str) -> OwnedRecord {
        OwnedRecord {
            head: head.to_string().as_bytes().to_vec(),
            seq: self
                .get_sequence()
                .into_iter()
                .map(|symbol| match symbol {
                    Some(s) => FASTA_ENCODE[&s],
                    None => 0x2d,
                })
                .collect(),
        }
    }

    pub fn get_length(&self) -> usize {
        match self {
            Haplotype::Wildtype(wt) => wt.get_length(),
            Haplotype::Descendant(ht) => ht.get_length(),
            Haplotype::Recombinant(rc) => rc.get_length(),
        }
    }

    pub fn get_tree(&self) -> String {
        let tree = self.get_subtree(self.get_reference().get_weak());
        format!("({});", tree)
    }

    pub fn get_subtree(&self, ancestor: HaplotypeWeak) -> String {
        self.prune_descendants();
        match self {
            Haplotype::Wildtype(wt) => wt.get_subtree(),
            Haplotype::Descendant(ht) => ht.get_subtree(),
            Haplotype::Recombinant(rc) => rc.get_subtree(ancestor),
        }
    }

    pub fn prune_tree(&self) {
        self.prune_descendants();
        self.get_descendants()
            .lock()
            .iter()
            .filter_map(|x| x.upgrade())
            .for_each(|x| x.prune_tree());
    }

    fn prune_descendants(&self) {
        let mut descendants_guard = self.get_descendants().lock();
        for idx in 0..descendants_guard.len() {
            if let Some(descendant) = descendants_guard[idx].upgrade() {
                let new = Haplotype::prune_descendant(descendant);
                descendants_guard[idx] = new.get_weak();
            }
        }
    }

    fn prune_descendant(descendant: HaplotypeRef) -> HaplotypeRef {
        let mut current = descendant;
        let mut chain: Vec<HaplotypeRef> = vec![current.clone()];

        // find chain of non-recombinant descendants
        loop {
            if !current.is_descendant() {
                break;
            }

            if current.n_descendants() != 1 {
                break;
            }

            let inner = current.unwrap_mutant();

            let next = inner
                .descendants
                .lock()
                .first()
                .unwrap()
                .clone()
                .upgrade()
                .unwrap();

            if !next.is_descendant() {
                break;
            }

            chain.push(next.clone());
            current = next;
        }

        // if there is nothing to prune, return
        if chain.len() == 1 {
            return current;
        }

        // aggregate changes
        let mut changes: HashMap<usize, (Symbol, Symbol)> = HashMap::new();
        for current in chain.iter().rev() {
            let inner = current.unwrap_mutant();
            inner
                .changes
                .iter()
                .for_each(|(position, change)| match changes.get_mut(position) {
                    Some((from, _to)) => *from = change.0,
                    None => {
                        changes.insert(*position, *change);
                    }
                });
        }

        // extract ancestor and wildtype
        let first = chain.first().unwrap();
        let ancestor = first.get_ancestors().0.unwrap();
        let wildtype = first.get_wildtype().get_weak();

        // extract descendants and generation
        let last = chain.last().unwrap();
        let generation = last.get_generation();

        // create new node
        unsafe {
            Descendant::new_and_replace(
                last.get_reference().clone(),
                ancestor.clone(),
                wildtype,
                changes,
                generation,
            );
        }

        last.clone()
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
                dirty_descendants: AtomicUsize::new(0),
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

impl Descendant {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        ancestor: HaplotypeRef,
        wildtype: HaplotypeWeak,
        changes: HashMap<usize, (Symbol, Symbol)>,
        generation: usize,
    ) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Descendant(Self {
                reference: reference.clone(),
                wildtype: wildtype.clone(),
                ancestor: ancestor.clone(),
                changes: changes.clone(),
                generation,
                fitness: OnceLock::new(),
                descendants: DescendantsCell::new(),
                dirty_descendants: AtomicUsize::new(0),
            })
        })
    }

    /// this is not thread-safe
    unsafe fn new_and_replace(
        old: HaplotypeRef,
        ancestor: HaplotypeRef,
        wildtype: HaplotypeWeak,
        changes: HashMap<usize, (Symbol, Symbol)>,
        generation: usize,
    ) {
        let descendants: Vec<HaplotypeWeak> = old
            .get_descendants()
            .lock()
            .iter()
            .cloned()
            .filter(|x| x.exists())
            .collect();
        let new = HaplotypeRef::new(Haplotype::Descendant(Self {
            reference: old.get_weak(),
            wildtype,
            ancestor,
            changes,
            generation,
            fitness: OnceLock::new(),
            descendants: DescendantsCell::from_iter(descendants),
            dirty_descendants: AtomicUsize::new(0),
        }));

        let old_ptr = old.as_ptr() as *mut Haplotype;
        let new_ptr = new.as_ptr() as *mut Haplotype;

        // evil hack to replace the old reference with the new one and keep
        // the reference count and reference pointer...
        // LOL
        std::ptr::swap(old_ptr, new_ptr);
    }

    pub fn get_base(&self, position: &usize) -> Symbol {
        match self.changes.get(position) {
            Some((_from, to)) => *to,
            None => self.ancestor.get_base(position),
        }
    }

    pub fn get_generation(&self) -> usize {
        self.generation
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    pub fn get_mutations(&self) -> HashMap<usize, (Symbol, Symbol)> {
        let mut mutations = self.ancestor.get_mutations();
        let wt_ref = self.wildtype.upgrade().unwrap();

        self.changes.iter().for_each(|(position, change)| {
            let wt_base = wt_ref.get_base(position);
            if change.1 == wt_base {
                mutations.remove(position);
            } else {
                mutations.insert(*position, *change);
            }
        });

        mutations
    }

    pub fn get_fitness(&self, fitness_table: &FitnessTable) -> &f64 {
        self.fitness.get_or_init(|| {
            let mut fitness = self.ancestor.get_fitness(fitness_table);

            self.changes.iter().for_each(|(position, change)| {
                fitness /= fitness_table.get_fitness(position, &change.0);
                fitness *= fitness_table.get_fitness(position, &change.1);
            });

            fitness_table.utility(fitness)
        })
    }

    pub fn get_subtree(&self) -> String {
        let current = self.reference.upgrade().unwrap();

        let name = current.get_string();
        let block_id = current.get_block_id();
        let branch_length = current.get_generation() - self.ancestor.get_generation();

        let node = format!("'{name}_{block_id}':{branch_length}");

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
                fitness: OnceLock::new(),
                descendants: DescendantsCell::new(),
                dirty_descendants: AtomicUsize::new(0),
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

    pub fn get_mutations(&self) -> HashMap<usize, (Symbol, Symbol)> {
        let left_mutations = self.left_ancestor.get_mutations();
        let right_mutations = self.right_ancestor.get_mutations();

        let mut mutations = HashMap::new();

        for (position, change) in left_mutations {
            if position >= self.left_position && position < self.right_position {
                mutations.insert(position, change);
            }
        }

        for (position, change) in right_mutations {
            if position < self.left_position || position >= self.right_position {
                mutations.insert(position, change);
            }
        }

        mutations
    }

    pub fn get_fitness(&self, fitness_table: &FitnessTable) -> &f64 {
        self.fitness.get_or_init(|| {
            let mutations = self.get_mutations();
            let mut fitness = 1.;

            for (position, (_from, to)) in mutations {
                fitness *= fitness_table.get_fitness(&position, &to);
            }

            fitness_table.utility(fitness)
        })
    }

    pub fn get_subtree(&self, ancestor: HaplotypeWeak) -> String {
        let name = self.reference.upgrade().unwrap().get_string();
        let block_id = self.reference.get_block_id();
        let branch_length = self.generation - ancestor.upgrade().unwrap().get_generation();

        let node = format!("#R'{name}_{block_id}':{branch_length}");

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
        assert_eq!(wt.get_tree(), format!("('0:0->3_{}':1)wt;", ht_id));

        let rc = Haplotype::create_recombinant(&wt, &ht, 1, 2, 2);
        let rc_id = rc.get_block_id();
        assert_eq!(
            wt.get_tree(),
            format!("((#R'0:0->3_{rc_id}':1)'0:0->3_{ht_id}':1,#R'0:0->3_{rc_id}':2)wt;")
        );
    }

    #[test]
    fn pruning_simple() {
        let bytes = vec![Some(0x00), Some(0x00), Some(0x00), Some(0x00)];
        let wt = Wildtype::new(bytes);

        let ht = wt.create_descendant(vec![0], vec![Some(0x01)], 1);
        let ht2 = ht.create_descendant(vec![0], vec![Some(0x02)], 2);
        let ht3 = ht2.create_descendant(vec![0], vec![Some(0x03)], 3);

        wt.prune_tree();

        assert_eq!(wt.n_descendants(), 1);

        let d = wt
            .get_descendants()
            .lock()
            .first()
            .unwrap()
            .upgrade()
            .unwrap();

        assert_eq!(d, ht3);
        assert_eq!(d.get_generation(), 3);
        assert_eq!(d.get_mutations().len(), 1);
        assert_eq!(d.get_mutations().get(&0), Some(&(Some(0x00), Some(0x03))));
    }
}
