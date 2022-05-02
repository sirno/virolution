use derivative::Derivative;
use std::cell::RefCell;
use std::ops::Range;
use std::rc::Rc;

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
    #[derivative(Debug(format_with = "print_reference_option"))]
    reference: Option<Rc<RefCell<Haplotype>>>,
    sequence: Vec<u8>,
    descendants: Vec<Rc<RefCell<Haplotype>>>,
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Descendant {
    #[derivative(Debug(format_with = "print_reference_option"))]
    reference: Option<Rc<RefCell<Haplotype>>>,
    #[derivative(Debug(format_with = "print_reference"))]
    wildtype: Rc<RefCell<Haplotype>>,
    #[derivative(Debug(format_with = "print_reference"))]
    ancestor: Rc<RefCell<Haplotype>>,
    descendants: Vec<Rc<RefCell<Haplotype>>>,
    position: usize,
    change: u8,
}

#[derive(Derivative)]
#[derivative(Debug)]
pub struct Recombinant {
    #[derivative(Debug(format_with = "print_reference_option"))]
    reference: Option<Rc<RefCell<Haplotype>>>,
    #[derivative(Debug(format_with = "print_reference"))]
    wildtype: Rc<RefCell<Haplotype>>,
    #[derivative(Debug(format_with = "print_reference"))]
    left_ancestor: Rc<RefCell<Haplotype>>,
    #[derivative(Debug(format_with = "print_reference"))]
    right_ancestor: Rc<RefCell<Haplotype>>,
    descendants: Vec<Rc<RefCell<Haplotype>>>,
    left_position: usize,
    right_position: usize,
}

fn print_reference(
    reference: &Rc<RefCell<Haplotype>>,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    write!(formatter, "{:x}", Rc::as_ptr(reference) as u64)
}

fn print_reference_option(
    reference_option: &Option<Rc<RefCell<Haplotype>>>,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    match reference_option {
        Some(reference) => write!(formatter, "{:x}", Rc::as_ptr(reference) as u64),
        None => write!(formatter, "None"),
    }
}

impl Haplotype {
    pub fn create_descendant(&mut self, position: usize, change: u8) -> Rc<RefCell<Haplotype>> {
        let ancestor = Rc::clone(self.get_reference());
        let wildtype = self.get_wildtype();

        let descendant = Rc::new(RefCell::new(Haplotype::Descendant(Descendant::new(
            ancestor, wildtype, position, change,
        ))));
        descendant
            .borrow_mut()
            .add_reference(Rc::clone(&descendant));

        self.add_descendant(&descendant);
        descendant
    }

    pub fn create_recombinant(
        left_ancestor: &Rc<RefCell<Haplotype>>,
        right_ancestor: &Rc<RefCell<Haplotype>>,
        left_position: usize,
        right_position: usize,
    ) -> Rc<RefCell<Haplotype>> {
        let wildtype = left_ancestor.borrow().get_wildtype();

        let recombinant = Rc::new(RefCell::new(Haplotype::Recombinant(Recombinant::new(
            wildtype,
            Rc::clone(left_ancestor),
            Rc::clone(right_ancestor),
            left_position,
            right_position,
        ))));
        recombinant
            .borrow_mut()
            .add_reference(Rc::clone(&recombinant));

        left_ancestor.borrow_mut().add_descendant(&recombinant);
        right_ancestor.borrow_mut().add_descendant(&recombinant);
        recombinant
    }

    pub fn add_reference(&mut self, reference: Rc<RefCell<Haplotype>>) {
        match self {
            Haplotype::Wildtype(wt) => wt.reference = Some(reference),
            Haplotype::Descendant(ht) => ht.reference = Some(reference),
            Haplotype::Recombinant(rc) => rc.reference = Some(reference),
        };
    }

    pub fn get_reference(&self) -> &Rc<RefCell<Haplotype>> {
        let option = match self {
            Haplotype::Wildtype(wt) => &wt.reference,
            Haplotype::Descendant(ht) => &ht.reference,
            Haplotype::Recombinant(rc) => &rc.reference,
        };
        match option {
            Some(reference) => &reference,
            None => {
                eprintln!("Haplotype incorrectly initialized.");
                std::process::exit(-1);
            }
        }
    }

    pub fn get_wildtype(&self) -> Rc<RefCell<Haplotype>> {
        match self {
            Haplotype::Wildtype(wt) => Rc::clone(&wt.get_reference()),
            Haplotype::Descendant(ht) => Rc::clone(&ht.wildtype),
            Haplotype::Recombinant(rc) => Rc::clone(&rc.wildtype),
        }
    }

    pub fn get_descendants(&self) -> &Vec<Rc<RefCell<Haplotype>>> {
        match self {
            Haplotype::Wildtype(wt) => &wt.descendants,
            Haplotype::Descendant(ht) => &ht.descendants,
            Haplotype::Recombinant(rc) => &rc.descendants,
        }
    }

    pub fn add_descendant(&mut self, descendant: &Rc<RefCell<Haplotype>>) {
        match self {
            Haplotype::Wildtype(wt) => wt.descendants.push(Rc::clone(descendant)),
            Haplotype::Descendant(ht) => ht.descendants.push(Rc::clone(descendant)),
            Haplotype::Recombinant(rc) => rc.descendants.push(Rc::clone(descendant)),
        };
    }

    pub fn get_base(&self, position: usize) -> u8 {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Descendant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_sequence(&self) -> Vec<u8> {
        let mut ranges: Vec<(Rc<RefCell<Haplotype>>, Range<usize>)> = Vec::new();
        ranges.push((Rc::clone(self.get_reference()), 0..self.get_length()));
        let mut sequence = vec![0; self.get_length()];
        while let Some((current, range)) = ranges.pop() {
            match &*(Rc::clone(&current)).borrow() {
                Haplotype::Wildtype(wt) => {
                    for (position, symbol) in wt.sequence.iter().enumerate() {
                        if range.contains(&position) && sequence[position] == 0 {
                            sequence[position] = *symbol;
                        }
                    }
                }
                Haplotype::Descendant(ht) => {
                    if range.contains(&ht.position) && sequence[ht.position] == 0 {
                        sequence[ht.position] = ht.change
                    }
                    ranges.push((Rc::clone(&ht.ancestor), range));
                }
                Haplotype::Recombinant(rc) => {
                    if range.end < rc.left_position || range.start > rc.right_position {
                        ranges.push((Rc::clone(&rc.right_ancestor), range))
                    } else if range.contains(&rc.left_position)
                        && range.contains(&rc.right_position)
                    {
                        ranges.push((Rc::clone(&rc.right_ancestor), range.start..rc.left_position));
                        ranges.push((
                            Rc::clone(&rc.left_ancestor),
                            rc.left_position..rc.right_position,
                        ));
                        ranges.push((Rc::clone(&rc.right_ancestor), rc.right_position..range.end));
                    } else if range.contains(&rc.left_position) {
                        ranges.push((Rc::clone(&rc.right_ancestor), range.start..rc.left_position));
                        ranges.push((Rc::clone(&rc.left_ancestor), rc.left_position..range.end));
                    } else if range.contains(&rc.right_position) {
                        ranges.push((Rc::clone(&rc.left_ancestor), range.start..rc.right_position));
                        ranges.push((Rc::clone(&rc.right_ancestor), rc.right_position..range.end));
                    } else {
                        ranges.push((Rc::clone(&rc.left_ancestor), range))
                    }
                }
            }
        }
        sequence
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
    pub fn create_wildtype(sequence: Vec<u8>) -> Rc<RefCell<Haplotype>> {
        let reference = Rc::new(RefCell::new(Haplotype::Wildtype(Self {
            reference: None,
            sequence: sequence,
            descendants: Vec::new(),
        })));
        reference.borrow_mut().add_reference(Rc::clone(&reference));
        reference
    }

    pub fn get_reference(&self) -> &Rc<RefCell<Haplotype>> {
        match &self.reference {
            Some(reference) => &reference,
            None => {
                eprintln!("Haplotype incorrectly initialized.");
                std::process::exit(-1);
            }
        }
    }

    pub fn get_base(&self, position: usize) -> u8 {
        self.sequence[position]
    }

    pub fn get_length(&self) -> usize {
        self.sequence.len()
    }
}

impl Descendant {
    pub fn new(
        ancestor: Rc<RefCell<Haplotype>>,
        wildtype: Rc<RefCell<Haplotype>>,
        position: usize,
        change: u8,
    ) -> Self {
        Self {
            reference: None,
            wildtype: wildtype,
            ancestor: ancestor,
            descendants: Vec::new(),
            position: position,
            change: change,
        }
    }

    pub fn get_reference(&self) -> &Rc<RefCell<Haplotype>> {
        match &self.reference {
            Some(reference) => &reference,
            None => {
                eprintln!("Haplotype incorrectly initialized.");
                std::process::exit(-1);
            }
        }
    }

    pub fn get_base(&self, position: usize) -> u8 {
        if self.position == position {
            return self.change;
        }
        self.ancestor.as_ref().borrow().get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.borrow().get_length()
    }
}

impl Recombinant {
    pub fn new(
        wildtype: Rc<RefCell<Haplotype>>,
        left_ancestor: Rc<RefCell<Haplotype>>,
        right_ancestor: Rc<RefCell<Haplotype>>,
        left_position: usize,
        right_position: usize,
    ) -> Self {
        Self {
            reference: None,
            wildtype: wildtype,
            left_ancestor: left_ancestor,
            right_ancestor: right_ancestor,
            left_position: left_position,
            right_position: right_position,
            descendants: Vec::new(),
        }
    }

    pub fn get_base(&self, position: usize) -> u8 {
        if position >= self.left_position && position < self.right_position {
            return self.left_ancestor.borrow().get_base(position);
        }

        self.right_ancestor.borrow().get_base(position)
    }

    pub fn get_length(&self) -> usize {
        self.wildtype.borrow().get_length()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn initiate_wildtype() {
        let bytes = vec![0x00, 0x01, 0x02, 0x03];
        let wt = Wildtype::create_wildtype(bytes);
        for i in 0..4 {
            assert_eq!(wt.borrow().get_base(i), i as u8);
        }
    }

    #[test]
    fn create_descendant() {
        let bytes = vec![0x00, 0x01, 0x02, 0x03];
        let wt = Wildtype::create_wildtype(bytes);
        let ht = wt.borrow_mut().create_descendant(0, 0x03);
        assert_eq!(ht.borrow().get_base(0), 0x03);
        assert_eq!(ht.borrow().get_base(1), 0x01);
        assert_eq!(ht.borrow().get_base(2), 0x02);
        assert_eq!(ht.borrow().get_base(3), 0x03);
    }

    #[test]
    fn create_wide_geneaology() {
        let bytes = vec![0x00, 0x01, 0x02, 0x03];
        let wt = Wildtype::create_wildtype(bytes);
        for i in 0..100 {
            let _ht = wt.borrow_mut().create_descendant(0, i);
        }
        for (position, descendant) in wt.borrow().get_descendants().iter().enumerate() {
            assert_eq!(descendant.borrow().get_base(0), position as u8)
        }
    }

    #[test]
    fn single_recombination() {
        let mut haplotypes = Vec::new();
        let bytes = vec![0x01; 100];
        let wildtype = Wildtype::create_wildtype(bytes);
        haplotypes.push(Rc::clone(&wildtype));
        for i in 0..100 {
            let ht = Rc::clone(haplotypes.last().unwrap());
            haplotypes.push(Rc::clone(&ht.borrow_mut().create_descendant(i, 0x02)));
        }
        haplotypes.push(Rc::clone(&wildtype.borrow_mut().create_descendant(0, 0x03)));
        for i in 1..100 {
            let ht = Rc::clone(haplotypes.last().unwrap());
            haplotypes.push(Rc::clone(&ht.borrow_mut().create_descendant(i, 0x03)));
        }
        let left_ancestor = &haplotypes[100];
        let right_ancestor = &haplotypes[200];
        let recombinant = Haplotype::create_recombinant(left_ancestor, right_ancestor, 10, 90);

        let mut expected = vec![0x03; 10];
        expected.append(&mut vec![0x02; 80]);
        expected.append(&mut vec![0x03; 10]);

        let out = recombinant.borrow().get_sequence();

        assert_eq!(out, expected)
    }
}
