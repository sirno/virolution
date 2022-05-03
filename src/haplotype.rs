use derivative::Derivative;
use std::cell::RefCell;
use std::fmt;
use std::ops::Range;
use std::rc::Rc;

type Symbol = Option<u8>;

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
    sequence: Vec<Symbol>,
    #[derivative(Debug(format_with = "print_descendants"))]
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
    #[derivative(Debug(format_with = "print_descendants"))]
    descendants: Vec<Rc<RefCell<Haplotype>>>,
    position: usize,
    change: Symbol,
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
    #[derivative(Debug(format_with = "print_descendants"))]
    descendants: Vec<Rc<RefCell<Haplotype>>>,
    left_position: usize,
    right_position: usize,
}

fn print_reference(
    reference: &Rc<RefCell<Haplotype>>,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    write!(formatter, "{}", reference.borrow().get_string())
}

fn print_reference_option(
    reference_option: &Option<Rc<RefCell<Haplotype>>>,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    match reference_option {
        Some(reference) => write!(formatter, "{}", reference.borrow().get_string()),
        None => write!(formatter, "None"),
    }
}

fn print_descendants(
    descendants: &Vec<Rc<RefCell<Haplotype>>>,
    formatter: &mut std::fmt::Formatter,
) -> Result<(), std::fmt::Error> {
    let out: Vec<String> = descendants
        .into_iter()
        .map(|descendant| descendant.borrow().get_string())
        .collect();
    write!(formatter, "{:?}", out)
}

impl fmt::Display for Haplotype {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.get_string())
    }
}

impl Haplotype {
    pub fn create_descendant(&mut self, position: usize, change: u8) -> Rc<RefCell<Haplotype>> {
        let ancestor = Rc::clone(self.get_reference());
        let wildtype = self.get_wildtype();

        let descendant = Rc::new(RefCell::new(Haplotype::Descendant(Descendant::new(
            ancestor,
            wildtype,
            position,
            Some(change),
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

    pub fn get_base(&self, position: usize) -> Symbol {
        match self {
            Haplotype::Wildtype(wt) => wt.get_base(position),
            Haplotype::Descendant(ht) => ht.get_base(position),
            Haplotype::Recombinant(rc) => rc.get_base(position),
        }
    }

    pub fn get_sequence(&self) -> Vec<Symbol> {
        let mut sequence = self.get_changes();
        if let Haplotype::Wildtype(wildtype) = &*Rc::clone(&self.get_wildtype()).borrow() {
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
        let mut ranges: Vec<(Rc<RefCell<Haplotype>>, Range<usize>)> = Vec::new();
        ranges.push((Rc::clone(self.get_reference()), 0..self.get_length()));
        let mut changes = vec![None; self.get_length()];
        while let Some((current, range)) = ranges.pop() {
            match &*(Rc::clone(&current)).borrow() {
                Haplotype::Wildtype(_wt) => {}
                Haplotype::Descendant(ht) => {
                    if let Some(symbol) = ht.change && range.contains(&ht.position) {
                        match changes[ht.position] {
                            Some(_) => {}
                            None => changes[ht.position] = Some(symbol),
                        }
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
        changes
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
    pub fn create_wildtype(sequence: Vec<Symbol>) -> Rc<RefCell<Haplotype>> {
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

    pub fn get_base(&self, position: usize) -> Symbol {
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
        change: Symbol,
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

    pub fn get_base(&self, position: usize) -> Symbol {
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

    pub fn get_base(&self, position: usize) -> Symbol {
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
        for i in 0..100 {
            let _ht = wt.borrow_mut().create_descendant(0, i);
        }
        for (position, descendant) in wt.borrow().get_descendants().iter().enumerate() {
            assert_eq!(descendant.borrow().get_base(0), Some(position as u8));
        }
    }

    #[test]
    fn single_recombination() {
        let mut haplotypes = Vec::new();
        let bytes = vec![Some(0x01); 100];
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

        let mut expected = vec![Some(0x03); 10];
        expected.append(&mut vec![Some(0x02); 80]);
        expected.append(&mut vec![Some(0x03); 10]);

        assert_eq!(recombinant.borrow().get_sequence(), expected);
    }

    #[test]
    fn get_length() {
        let wildtype = Wildtype::create_wildtype(vec![Some(0x00); 100]);
        let haplotype = wildtype.borrow_mut().create_descendant(0, 0x01);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75);
        assert_eq!(wildtype.borrow().get_length(), 100);
        assert_eq!(haplotype.borrow().get_length(), 100);
        assert_eq!(recombinant.borrow().get_length(), 100);
    }

    #[test]
    fn get_sequence() {
        let wildtype = Wildtype::create_wildtype(vec![Some(0x00); 100]);
        let haplotype = wildtype.borrow_mut().create_descendant(0, 0x01);
        let recombinant = Haplotype::create_recombinant(&wildtype, &haplotype, 25, 75);

        let mut expected = vec![Some(0x01)];
        expected.append(&mut vec![Some(0x00); 99]);

        assert_eq!(recombinant.borrow().get_sequence(), expected);
    }
}
