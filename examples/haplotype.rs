extern crate virolution;

use virolution::core::haplotype::*;

fn main() {
    let bytes = vec![Some(0x00); 4];
    let wt = Wildtype::new(bytes);
    let ht = wt.create_descendant(vec![2], vec![Some(0x01)], 0);
    let ht2 = ht.create_descendant(vec![1], vec![Some(0x02)], 0);
    let ht3 = ht2.create_descendant(vec![2], vec![Some(0x03)], 0);
    let ht4 = Haplotype::create_recombinant(&ht, &ht3, 0, 2, 0);
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
            .map(|idx| wt.get_base(&idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht.get_base(&idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht2.get_base(&idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht3.get_base(&idx))
            .collect::<Vec<Symbol>>()
    );
    println!(
        "{:?}",
        (0..4)
            .into_iter()
            .map(|idx| ht4.get_base(&idx))
            .collect::<Vec<Symbol>>()
    );

    println!("---get_sequence---");
    println!("{:?}", wt.get_sequence());
    println!("{:?}", ht.get_sequence());
    println!("{:?}", ht2.get_sequence());
    println!("{:?}", ht3.get_sequence());
    println!("{:?}", ht4.get_sequence());
}
