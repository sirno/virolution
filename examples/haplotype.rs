extern crate virolution;

use virolution::haplotype::*;

fn main() {
    let bytes = vec![Some(0x00); 4];
    let wt = Wildtype::create_wildtype(bytes);
    let ht = wt.borrow_mut().create_descendant(2, 0x01);
    let ht2 = ht.borrow_mut().create_descendant(1, 0x02);
    let ht3 = ht2.borrow_mut().create_descendant(2, 0x03);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2);
    println!("---debug---");
    println!("wt: {:?}", *wt);
    println!("ht: {:?}", *ht);
    println!("ht2: {:?}", *ht2);
    println!("ht3: {:?}", *ht3);
    println!("ht4: {:?}", *ht4);
    println!("---get_base---");
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| wt.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht2.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht3.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht4.borrow().get_base(idx))
            .collect::<Vec<Symbol>>()
    );

    println!("---get_sequence---");
    println!("{:?}", wt.borrow().get_sequence());
    println!("{:?}", ht.borrow().get_sequence());
    println!("{:?}", ht2.borrow().get_sequence());
    println!("{:?}", ht3.borrow().get_sequence());
    println!("{:?}", ht4.borrow().get_sequence());
}
