use super::fitness::FitnessTable;
use super::references::sync::{HaplotypeRef, HaplotypeWeak};
use derivative::Derivative;
use phf::phf_map;
use seq_io::fasta::OwnedRecord;
use std::collections::HashSet;
use std::fmt;
use std::ops::Range;

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

#[derive(Derivative)]
#[derivative(Debug)]
pub enum Haplotype {
    Wildtype(Wildtype),
    Descendant(Descendant),
    Recombinant(Recombinant),
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Wildtype {
    #[derivative(Debug(format_with = "print_reference_weak"))]
    reference: HaplotypeWeak,
    #[derivative(Debug = "ignore")]
    sequence: Vec<Symbol>,
    #[derivative(Debug(format_with = "print_descendants"))]
    descendants: Vec<HaplotypeWeak>,
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Descendant {
    #[derivative(Debug(format_with = "print_reference_weak"))]
    reference: HaplotypeWeak,
    #[derivative(Debug(format_with = "print_reference"))]
    wildtype: HaplotypeRef,
    #[derivative(Debug(format_with = "print_reference"))]
    ancestor: HaplotypeRef,
    #[derivative(Debug(format_with = "print_descendants"))]
    descendants: Vec<HaplotypeWeak>,
    position: usize,
    change: Symbol,
    fitness: Option<f64>,
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Recombinant {
    #[derivative(Debug(format_with = "print_reference_weak"))]
    reference: HaplotypeWeak,
    #[derivative(Debug(format_with = "print_reference"))]
    wildtype: HaplotypeRef,
    #[derivative(Debug(format_with = "print_reference"))]
    left_ancestor: HaplotypeRef,
    #[derivative(Debug(format_with = "print_reference"))]
    right_ancestor: HaplotypeRef,
    #[derivative(Debug(format_with = "print_descendants"))]
    descendants: Vec<HaplotypeWeak>,
    left_position: usize,
    right_position: usize,
    fitness: Option<f64>,
}

fn print_reference(
    reference: &HaplotypeRef,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    write!(formatter, "{}", reference.borrow().get_string())
}

fn print_reference_weak(
    reference_weak: &HaplotypeWeak,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    match reference_weak.upgrade() {
        Some(reference) => write!(formatter, "{}", reference.borrow().get_string()),
        None => write!(formatter, "None"),
    }
}

fn print_descendants(
    descendants: &Vec<HaplotypeWeak>,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    let out: Vec<String> = descendants
        .into_iter()
        .map(|descendant_weak| match descendant_weak.upgrade() {
            Some(descendant) => descendant.borrow().to_string(),
            None => "None".to_string(),
        })
        .collect();
    write!(formatter, "{:?}", out)
}

impl Haplotype {
    pub fn create_descendant(&mut self, position: usize, change: u8) -> HaplotypeRef {
        let ancestor = self.get_reference();
        let wildtype = self.get_wildtype();
        let descendant = Descendant::new(ancestor, wildtype, position, Some(change));

        self.add_descendant(descendant.get_weak());

        descendant
    }

    pub fn create_recombinant(
        left_ancestor: &HaplotypeRef,
        right_ancestor: &HaplotypeRef,
        left_position: usize,
        right_position: usize,
    ) -> HaplotypeRef {
        let wildtype = left_ancestor.borrow().get_wildtype();

        let recombinant = Recombinant::new(
            wildtype,
            left_ancestor.get_clone(),
            right_ancestor.get_clone(),
            left_position,
            right_position,
        );

        left_ancestor
            .borrow_mut()
            .add_descendant(recombinant.get_weak());
        right_ancestor
            .borrow_mut()
            .add_descendant(recombinant.get_weak());
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
            Haplotype::Descendant(ht) => ht.wildtype.get_clone(),
            Haplotype::Recombinant(rc) => rc.wildtype.get_clone(),
        }
    }

    pub fn get_descendants(&self) -> &Vec<HaplotypeWeak> {
        match self {
            Haplotype::Wildtype(wt) => &wt.descendants,
            Haplotype::Descendant(ht) => &ht.descendants,
            Haplotype::Recombinant(rc) => &rc.descendants,
        }
    }

    fn add_descendant(&mut self, descendant: HaplotypeWeak) {
        match self {
            Haplotype::Wildtype(wt) => wt.descendants.push(descendant),
            Haplotype::Descendant(ht) => ht.descendants.push(descendant),
            Haplotype::Recombinant(rc) => rc.descendants.push(descendant),
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
        let mut sequence = self.get_changes();
        if let Haplotype::Wildtype(wildtype) = &*self.get_wildtype().get_clone().borrow() {
            for (position, symbol) in wildtype.sequence.iter().enumerate() {
                match sequence[position] {
                    None => sequence[position] = *symbol,
                    Some(_) => {}
                }
            }
        }
        sequence
    }

    pub fn get_wildtype_sequence(&self) -> Vec<Symbol> {
        match self {
            Haplotype::Wildtype(wt) => wt.sequence.to_vec(),
            Haplotype::Descendant(_ht) => self.get_wildtype().borrow().get_wildtype_sequence(),
            Haplotype::Recombinant(_rc) => self.get_wildtype().borrow().get_wildtype_sequence(),
        }
    }

    pub fn get_string(&self) -> String {
        let changes = self.get_changes();
        let wildtype = self.get_wildtype_sequence();

        let mut out = String::new();
        for (position, change) in changes.iter().enumerate() {
            match change {
                Some(symbol) => {
                    if let Some(wt_symbol) = wildtype[position] {
                        out.push_str(format!(";{position}:{wt_symbol}->{symbol}").as_str())
                    }
                }
                None => {}
            }
        }

        if out.len() == 0 {
            return "wt".to_string();
        }

        out.remove(0);
        out
    }

    pub fn get_changes(&self) -> Vec<Symbol> {
        let mut ranges: Vec<(HaplotypeRef, Range<usize>)> = Vec::new();
        ranges.push((self.get_reference().get_clone(), 0..self.get_length()));
        let mut changes = vec![None; self.get_length()];
        while let Some((current, range)) = ranges.pop() {
            match &*current.get_clone().borrow() {
                Haplotype::Wildtype(_wt) => {}
                Haplotype::Descendant(ht) => {
                    if let Some(symbol) = ht.change && range.contains(&ht.position) {
                        match changes[ht.position] {
                            Some(_) => {},
                            None => changes[ht.position] = Some(symbol),
                        }
                    }
                    ranges.push((ht.ancestor.get_clone(), range));
                }
                Haplotype::Recombinant(rc) => rc.push_to_ranges(&mut ranges, range),
            }
        }
        changes
    }

    pub fn get_fitness(&self, fitness_table: &FitnessTable) -> f64 {
        match self {
            Haplotype::Wildtype(_wt) => return 1.,
            Haplotype::Descendant(Descendant {
                fitness: Some(val), ..
            }) => return *val,
            Haplotype::Recombinant(Recombinant {
                fitness: Some(val), ..
            }) => return *val,
            _ => {}
        }

        let mut ranges: Vec<(HaplotypeRef, Range<usize>)> = Vec::new();
        let mut checked: HashSet<usize> = HashSet::new();
        ranges.push((self.get_reference().get_clone(), 0..self.get_length()));
        let mut fitness = 1.;
        while let Some((current, range)) = ranges.pop() {
            match &*current.get_clone().borrow() {
                Haplotype::Wildtype(_wt) => {}
                Haplotype::Descendant(ht) => {
                    if range.contains(&ht.position) && !checked.contains(&ht.position) {
                        fitness *= fitness_table.get_fitness(&ht.position, &ht.change);
                        checked.insert(ht.position);
                    }
                    ranges.push((ht.ancestor.get_clone(), range));
                }
                Haplotype::Recombinant(rc) => rc.push_to_ranges(&mut ranges, range),
            }
        }
        fitness
    }

    pub fn get_record(&self, head: &str) -> OwnedRecord {
        let signature = self.to_string();
        let header = format!("{head};haplotype={signature}");
        OwnedRecord {
            head: header.as_bytes().to_vec(),
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
}

impl Wildtype {
    pub fn create_wildtype(sequence: Vec<Symbol>) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Wildtype(Self {
                reference: reference.clone(),
                sequence: sequence.clone(),
                descendants: Vec::new(),
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
}

impl Descendant {
    pub fn new(
        ancestor: HaplotypeRef,
        wildtype: HaplotypeRef,
        position: usize,
        change: Symbol,
    ) -> HaplotypeRef {
        HaplotypeRef::new_cyclic(|reference| {
            Haplotype::Descendant(Self {
                reference: reference.clone(),
                wildtype: wildtype.get_clone(),
                ancestor: ancestor.get_clone(),
                descendants: Vec::new(),
                position: position,
                change: change,
                fitness: None,
            })
        })
    }

    pub fn get_base(&self, position: usize) -> Symbol {
        if self.position == position {
            return self.change;
        }
        self.ancestor.borrow().get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.borrow().get_length()
    }
}

impl Recombinant {
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
                wildtype: wildtype.get_clone(),
                left_ancestor: left_ancestor.get_clone(),
                right_ancestor: right_ancestor.get_clone(),
                left_position: left_position,
                right_position: right_position,
                descendants: Vec::new(),
                fitness: None,
            })
        })
    }

    pub fn get_base(&self, position: usize) -> Symbol {
        if position >= self.left_position && position < self.right_position {
            return self.left_ancestor.borrow().get_base(position);
        }

        self.right_ancestor.borrow().get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.borrow().get_length()
    }

    fn push_to_ranges(&self, ranges: &mut Vec<(HaplotypeRef, Range<usize>)>, range: Range<usize>) {
        if range.end < self.left_position || range.start > self.right_position {
            ranges.push((self.right_ancestor.get_clone(), range))
        } else if range.contains(&self.left_position) && range.contains(&self.right_position) {
            ranges.push((
                self.right_ancestor.get_clone(),
                range.start..self.left_position,
            ));
            ranges.push((
                self.left_ancestor.get_clone(),
                self.left_position..self.right_position,
            ));
            ranges.push((
                self.right_ancestor.get_clone(),
                self.right_position..range.end,
            ));
        } else if range.contains(&self.left_position) {
            ranges.push((
                self.right_ancestor.get_clone(),
                range.start..self.left_position,
            ));
            ranges.push((
                self.left_ancestor.get_clone(),
                self.left_position..range.end,
            ));
        } else if range.contains(&self.right_position) {
            ranges.push((
                self.left_ancestor.get_clone(),
                range.start..self.right_position,
            ));
            ranges.push((
                self.right_ancestor.get_clone(),
                self.right_position..range.end,
            ));
        } else {
            ranges.push((self.left_ancestor.get_clone(), range))
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

    #[test]
    fn initiate_wildtype() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::create_wildtype(bytes);
        for i in 0..4 {
            assert_eq!(wt.borrow().get_base(i), Some(i as u8));
        }
    }

    #[test]
    fn create_descendant() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::create_wildtype(bytes);
        let ht = wt.borrow_mut().create_descendant(0, 0x03);
        assert_eq!(ht.borrow().get_base(0), Some(0x03));
        assert_eq!(ht.borrow().get_base(1), Some(0x01));
        assert_eq!(ht.borrow().get_base(2), Some(0x02));
        assert_eq!(ht.borrow().get_base(3), Some(0x03));
    }

    #[test]
    fn create_wide_geneaology() {
        let bytes = vec![Some(0x00), Some(0x01), Some(0x02), Some(0x03)];
        let wt = Wildtype::create_wildtype(bytes);
        let _hts: Vec<HaplotypeRef> = (0..100)
            .map(|i| wt.borrow_mut().create_descendant(0, i))
            .collect();
        for (position, descendant) in wt.borrow().get_descendants().iter().enumerate() {
            if let Some(d) = descendant.upgrade() {
                assert_eq!(d.borrow().get_base(0), Some(position as u8));
            } else {
                panic!();
            }
        }
    }

    #[test]
    fn single_recombination() {
        let mut haplotypes: Vec<HaplotypeRef> = Vec::new();
        let bytes = vec![Some(0x01); 100];
        let wildtype = Wildtype::create_wildtype(bytes);
        haplotypes.push(wildtype.get_clone());
        for i in 0..100 {
            let ht = haplotypes.last().unwrap().get_clone();
            haplotypes.push(ht.borrow_mut().create_descendant(i, 0x02));
        }
        haplotypes.push(wildtype.borrow_mut().create_descendant(0, 0x03));
        for i in 1..100 {
            let ht = haplotypes.last().unwrap().get_clone();
            haplotypes.push(ht.borrow_mut().create_descendant(i, 0x03));
        }
        let mut left_ancestor = haplotypes[100].get_clone();
        let mut right_ancestor = haplotypes[200].get_clone();
        let recombinant =
            Haplotype::create_recombinant(&mut left_ancestor, &mut right_ancestor, 10, 90);

        let mut expected = vec![Some(0x03); 10];
        expected.append(&mut vec![Some(0x02); 80]);
        expected.append(&mut vec![Some(0x03); 10]);

        assert_eq!(recombinant.borrow().get_sequence(), expected);
    }

    #[test]
    fn get_length() {
        let mut wildtype = Wildtype::create_wildtype(vec![Some(0x00); 100]);
        let mut haplotype = wildtype.borrow_mut().create_descendant(0, 0x01);
        let recombinant = Haplotype::create_recombinant(&mut wildtype, &mut haplotype, 25, 75);
        assert_eq!(wildtype.borrow().get_length(), 100);
        assert_eq!(haplotype.borrow().get_length(), 100);
        assert_eq!(recombinant.borrow().get_length(), 100);
    }

    #[test]
    fn get_sequence() {
        let mut wildtype = Wildtype::create_wildtype(vec![Some(0x00); 100]);
        let mut haplotype = wildtype.borrow_mut().create_descendant(0, 0x01);
        let recombinant = Haplotype::create_recombinant(&mut wildtype, &mut haplotype, 25, 75);

        let mut expected = vec![Some(0x01)];
        expected.append(&mut vec![Some(0x00); 99]);

        assert_eq!(recombinant.borrow().get_sequence(), expected);
    }
}
