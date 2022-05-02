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

        match self {
            Haplotype::Wildtype(wt) => wt.descendants.push(Rc::clone(&descendant)),
            Haplotype::Descendant(ht) => ht.descendants.push(Rc::clone(&descendant)),
            Haplotype::Recombinant(rc) => rc.descendants.push(Rc::clone(&descendant)),
        };
        descendant
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

    pub fn get_base(&self, position: usize) -> u8 {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Descendant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_sequence(&self) -> Vec<u8> {
        // let mut ranges: Vec<(Rc<RefCell<Haplotype>>, Range<usize>)> = Vec::new();
        let mut sequence = vec![0; self.get_length()];
        let mut current = Rc::clone(self.get_reference());
        while let Haplotype::Descendant(ht) = &*(Rc::clone(&current)).borrow() {
            if sequence[ht.position] == 0 {
                sequence[ht.position] = ht.change
            }
            current = Rc::clone(&ht.ancestor);
        }
        if let Haplotype::Wildtype(wt) = &*(*current).borrow() {
            for (position, symbol) in wt.sequence.iter().enumerate() {
                if sequence[position] == 0 {
                    sequence[position] = *symbol;
                }
            }
        }
        return sequence;
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
            ancestor: ancestor,
            wildtype: wildtype,
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
}
