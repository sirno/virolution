use super::fitness::FitnessTable;
use super::references::sync::{HaplotypeRef, HaplotypeWeak};
use derivative::Derivative;
use phf::phf_map;
use seq_io::fasta::OwnedRecord;
use std::collections::HashMap;
use std::fmt;
use std::sync::{Arc, Mutex, OnceLock};

// #[derive(Clone, Debug, Deref)]
// pub struct Symbol(Option<u8>);
pub type Symbol = Option<u8>;

pub static FASTA_ENCODE: phf::Map<u8, u8> = phf_map! {
    0x00u8 => 0x41,
    0x01u8 => 0x43,
    0x02u8 => 0x47,
    0x03u8 => 0x54,
};
pub static FASTA_DECODE: phf::Map<u8, u8> = phf_map! {
    0x41u8 => 0x00,
    0x43u8 => 0x01,
    0x47u8 => 0x02,
    0x54u8 => 0x03,
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
    descendants: Arc<Mutex<Vec<HaplotypeWeak>>>,
}

#[derive(Debug)]
pub struct Descendant {
    reference: HaplotypeWeak,
    wildtype: HaplotypeWeak,
    ancestor: HaplotypeRef,
    descendants: Arc<Mutex<Vec<HaplotypeWeak>>>,
    position: usize,
    change: (Symbol, Symbol),
    changes: OnceLock<HashMap<usize, (Symbol, Symbol)>>,
    fitness: OnceLock<f64>,
}

#[derive(Debug)]
pub struct Recombinant {
    reference: HaplotypeWeak,
    wildtype: HaplotypeWeak,
    left_ancestor: HaplotypeRef,
    right_ancestor: HaplotypeRef,
    descendants: Arc<Mutex<Vec<HaplotypeWeak>>>,
    left_position: usize,
    right_position: usize,
    changes: OnceLock<HashMap<usize, (Symbol, Symbol)>>,
    fitness: OnceLock<f64>,
}

impl Haplotype {
    pub fn create_descendant(&self, position: usize, change: u8) -> HaplotypeRef {
        let ancestor = self.get_reference();
        let wildtype = self.get_wildtype();

        let descendant = Descendant::new(
            ancestor,
            wildtype,
            position,
            (self.get_base(position), Some(change)),
        );

        self.add_descendant(descendant.get_weak());

        descendant
    }

    pub fn create_recombinant(
        left_ancestor: &HaplotypeRef,
        right_ancestor: &HaplotypeRef,
        left_position: usize,
        right_position: usize,
    ) -> HaplotypeRef {
        let wildtype = left_ancestor.get_wildtype();

        let recombinant = Recombinant::new(
            wildtype,
            left_ancestor.get_clone(),
            right_ancestor.get_clone(),
            left_position,
            right_position,
        );

        left_ancestor.add_descendant(recombinant.get_weak());
        right_ancestor.add_descendant(recombinant.get_weak());
        recombinant
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
            // Haplotype::Wildtype(wt) => wt.get_reference().get_clone(),
            Haplotype::Wildtype(wt) => wt.get_reference().get_clone(),
            Haplotype::Descendant(ht) => ht.wildtype.upgrade().unwrap().get_clone(),
            Haplotype::Recombinant(rc) => rc.wildtype.upgrade().unwrap().get_clone(),
        }
    }

    pub fn get_descendants(&self) -> Arc<Mutex<Vec<HaplotypeWeak>>> {
        match self {
            Haplotype::Wildtype(wt) => wt.descendants.clone(),
            Haplotype::Descendant(ht) => ht.descendants.clone(),
            Haplotype::Recombinant(rc) => rc.descendants.clone(),
        }
    }

    fn add_descendant(&self, descendant: HaplotypeWeak) {
        match self {
            Haplotype::Wildtype(wt) => wt.descendants.lock().unwrap().push(descendant),
            Haplotype::Descendant(ht) => ht.descendants.lock().unwrap().push(descendant),
            Haplotype::Recombinant(rc) => rc.descendants.lock().unwrap().push(descendant),
        };
    }

    pub fn get_base(&self, position: usize) -> Symbol {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Descendant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_sequence(&self) -> Vec<Symbol> {
        let changes = self.get_changes();
        let mut sequence = self.get_wildtype_sequence();
        for (position, (_, to)) in changes {
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
        let changes = self.get_changes();
        let wildtype = self.get_wildtype_sequence();

        let mut out = String::new();
        for (position, (from, to)) in changes.iter() {
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

    pub fn get_changes(&self) -> HashMap<usize, (Symbol, Symbol)> {
        match self {
            Haplotype::Wildtype(_wt) => HashMap::new(),
            Haplotype::Descendant(ht) => ht.get_changes().clone(),
            Haplotype::Recombinant(rc) => rc.get_changes().clone(),
        }
    }

    pub fn get_fitness(&self, fitness_table: &FitnessTable) -> f64 {
        let compute_fitness = || {
            let changes = self.get_changes();
            let mut fitness = 1.;

            for (position, (_from, to)) in changes {
                fitness *= fitness_table.get_fitness(&position, &to);
            }

            fitness
        };

        match self {
            Haplotype::Wildtype(_wt) => 1.,
            Haplotype::Descendant(Descendant { fitness, .. }) => {
                *fitness.get_or_init(compute_fitness)
            }
            Haplotype::Recombinant(Recombinant { fitness, .. }) => {
                *fitness.get_or_init(compute_fitness)
            }
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
        let mut tree = match self {
            Haplotype::Wildtype(wt) => wt.get_subtree(),
            Haplotype::Descendant(ht) => ht.get_subtree(),
            Haplotype::Recombinant(rc) => rc.get_subtree(self.get_reference().get_weak()),
        };
        tree.push(';');
        tree
    }

    pub fn get_subtree(&self, ancestor: HaplotypeWeak) -> String {
        match self {
            Haplotype::Wildtype(wt) => wt.get_subtree(),
            Haplotype::Descendant(ht) => ht.get_subtree(),
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
                descendants: Arc::new(Mutex::new(Vec::new())),
            })
        })
    }

    pub fn get_reference(&self) -> HaplotypeRef {
        self.reference
            .upgrade()
            .expect("Self-reference has been dropped.")
    }

    pub fn get_base(&self, position: usize) -> Symbol {
        self.sequence[position]
    }

    pub fn get_length(&self) -> usize {
        self.sequence.len()
    }

    #[inline]
    pub fn get_subtree(&self) -> String {
        let inner = self
            .descendants
            .lock()
            .unwrap()
            .iter()
            .filter(|x| x.exists())
            .map(|x| x.upgrade().unwrap().get_subtree(self.reference.clone()))
            .collect::<Vec<String>>()
            .join(",");
        format!("({inner})wt")
    }
}

impl Descendant {
    #[allow(clippy::new_ret_no_self)]
    pub fn new(
        ancestor: HaplotypeRef,
        wildtype: HaplotypeRef,
        position: usize,
        change: (Symbol, Symbol),
    ) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Descendant(Self {
                reference: reference.clone(),
                wildtype: wildtype.get_weak(),
                ancestor: ancestor.get_clone(),
                descendants: Arc::new(Mutex::new(Vec::new())),
                position,
                change,
                changes: OnceLock::new(),
                fitness: OnceLock::new(),
            })
        })
    }

    pub fn get_base(&self, position: usize) -> Symbol {
        if self.position == position {
            return self.change.1;
        }
        self.ancestor.get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    pub fn get_changes(&self) -> &HashMap<usize, (Symbol, Symbol)> {
        self.changes.get_or_init(|| {
            let mut changes = self.ancestor.get_changes();
            let wt_base = self.wildtype.upgrade().unwrap().get_base(self.position);

            if self.change.1 == wt_base {
                changes.remove(&self.position);
                return changes;
            }

            changes.insert(self.position, (wt_base, self.change.1));
            changes
        })
    }

    #[inline]
    pub fn get_subtree(&self) -> String {
        let inner = self
            .descendants
            .lock()
            .unwrap()
            .iter()
            .filter(|x| x.exists())
            .map(|x| x.upgrade().unwrap().get_subtree(self.reference.clone()))
            .collect::<Vec<String>>()
            .join(",");
        let outer = self.reference.upgrade().unwrap().get_string();

        if inner.is_empty() {
            format!("'{outer}m{}'", self.reference.get_id())
        } else {
            format!("({inner})'{outer}m{}'", self.reference.get_id())
        }
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
    ) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Recombinant(Self {
                reference: reference.clone(),
                wildtype: wildtype.get_weak(),
                left_ancestor: left_ancestor.get_clone(),
                right_ancestor: right_ancestor.get_clone(),
                left_position,
                right_position,
                descendants: Arc::new(Mutex::new(Vec::new())),
                changes: OnceLock::new(),
                fitness: OnceLock::new(),
            })
        })
    }

    pub fn get_base(&self, position: usize) -> Symbol {
        if position >= self.left_position && position < self.right_position {
            return self.left_ancestor.get_base(position);
        }

        self.right_ancestor.get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.upgrade().unwrap().get_length()
    }

    pub fn get_changes(&self) -> &HashMap<usize, (Symbol, Symbol)> {
        self.changes.get_or_init(|| {
            let left_changes = self.left_ancestor.get_changes();
            let right_changes = self.right_ancestor.get_changes();

            let mut changes = HashMap::new();

            for (position, change) in left_changes {
                if position >= self.left_position && position < self.right_position {
                    changes.insert(position, change);
                }
            }

            for (position, change) in right_changes {
                if position < self.left_position || position >= self.right_position {
                    changes.insert(position, change);
                }
            }

            changes
        })
    }

    #[inline]
    pub fn get_subtree(&self, ancestor: HaplotypeWeak) -> String {
        let inner = if self.left_ancestor.get_weak() == ancestor {
            self.descendants
                .lock()
                .unwrap()
                .iter()
                .filter(|x| x.exists())
                .map(|x| x.upgrade().unwrap().get_subtree(self.reference.clone()))
                .collect::<Vec<String>>()
                .join(",")
        } else {
            "".to_string()
        };
        if inner.is_empty() {
            format!(
                "#R'{}r{}'",
                self.reference.upgrade().unwrap().get_string(),
                self.reference.get_id(),
            )
        } else {
            format!(
                "({})#R'{}r{}'",
                inner,
                self.reference.upgrade().unwrap().get_string(),
                self.reference.get_id(),
            )
        }
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
            assert_eq!(wt.get_base(i), Some(i as u8));
        }
    }

    #[test]
    fn create_descendant() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);
        let ht = wt.create_descendant(0, 0x03);
        assert_eq!(ht.get_base(0), Some(0x03));
        assert_eq!(ht.get_base(1), Some(0x01));
        assert_eq!(ht.get_base(2), Some(0x02));
        assert_eq!(ht.get_base(3), Some(0x03));
    }

    #[test]
    fn create_wide_geneaology() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);
        let _hts: Vec<HaplotypeRef> = (0..100).map(|i| wt.create_descendant(0, i)).collect();
        for (position, descendant) in wt.get_descendants().lock().unwrap().iter().enumerate() {
            if let Some(d) = descendant.upgrade() {
                assert_eq!(d.get_base(0), Some(position as u8));
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
        haplotypes.push(wildtype.get_clone());
        for i in 0..100 {
            let ht = haplotypes.last().unwrap().get_clone();
            haplotypes.push(ht.create_descendant(i, 0x02));
        }
        haplotypes.push(wildtype.create_descendant(0, 0x03));
        for i in 1..100 {
            let ht = haplotypes.last().unwrap().get_clone();
            haplotypes.push(ht.create_descendant(i, 0x03));
        }
        let left_ancestor = haplotypes[100].get_clone();
        let right_ancestor = haplotypes[200].get_clone();
        let recombinant = Haplotype::create_recombinant(&left_ancestor, &right_ancestor, 10, 90);

        let mut expected = vec![Some(0x03); 10];
        expected.append(&mut vec![Some(0x02); 80]);
        expected.append(&mut vec![Some(0x03); 10]);

        assert_eq!(recombinant.get_sequence(), expected);
    }

    #[test]
    #[serial]
    fn get_length() {
        let wildtype = Wildtype::new(vec![Some(0x00); 100]);
        let haplotype = wildtype.create_descendant(0, 0x01);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75);
        assert_eq!(wildtype.get_length(), 100);
        assert_eq!(haplotype.get_length(), 100);
        assert_eq!(recombinant.get_length(), 100);
    }

    #[test]
    #[serial]
    fn get_sequence() {
        let wildtype = Wildtype::new(vec![Some(0x00); 100]);
        let haplotype = wildtype.create_descendant(0, 0x01);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75);

        let mut expected = vec![Some(0x01)];
        expected.append(&mut vec![Some(0x00); 99]);

        assert_eq!(recombinant.get_sequence(), expected);
    }

    #[test]
    #[serial]
    fn create_tree() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::new(bytes);

        let ht = wt.create_descendant(0, 0x03);
        let ht_id = ht.get_id();
        assert_eq!(wt.get_tree(), format!("('0:0->3m{}')wt;", ht_id));

        let rc = Haplotype::create_recombinant(&wt, &ht, 1, 2);
        let rc_id = rc.get_id();
        assert_eq!(
            wt.get_tree(),
            format!("((#R'0:0->3r{rc_id}')'0:0->3m{ht_id}',#R'0:0->3r{rc_id}')wt;")
        );
    }
}
