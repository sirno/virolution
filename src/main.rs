#![feature(let_chains)]

mod haplotype;

use haplotype::*;

fn _haplotype_experiments() {
    let bytes = vec![Some(0x00); 4];
    let wt = Wildtype::create_wildtype(bytes);
    let ht = wt.borrow_mut().create_descendant(2, 0x01);
    let ht2 = ht.borrow_mut().create_descendant(1, 0x02);
    let ht3 = ht2.borrow_mut().create_descendant(2, 0x03);
    println!("wt: {:?}", *wt.borrow());
    println!("ht: {:?}", *ht.borrow());
    println!("ht2: {:?}", *ht2.borrow());
    println!("ht3: {:?}", *ht3.borrow());
    for i in 0..4 {
        println!("{}: {:?}", i, ht3.borrow().get_base(i));
    }
    let seq = ht2.borrow().get_sequence();
    println!("{:?}", seq);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2);
    println!("{:?}", ht4.borrow().get_sequence());
}

fn _population_experiments() {}

fn main() {
    _haplotype_experiments();
    _population_experiments();
}
