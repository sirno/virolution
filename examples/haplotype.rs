extern crate virolution;

use virolution::core::attributes::AttributeSetDefinition;
use virolution::core::haplotype::*;
use virolution::encoding::Nucleotide as Nt;

fn main() {
    let bytes = vec![Nt::A; 4];
    let attribute_definition = AttributeSetDefinition::new();
    let wt = Wildtype::new(bytes, &attribute_definition);
    let ht = wt.create_descendant(vec![2], vec![Nt::T], 0);
    let ht2 = ht.create_descendant(vec![1], vec![Nt::C], 0);
    let ht3 = ht2.create_descendant(vec![2], vec![Nt::G], 0);
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
        (0..4).map(|idx| wt.get_base(&idx)).collect::<Vec<Nt>>()
    );
    println!(
        "{:?}",
        (0..4).map(|idx| ht.get_base(&idx)).collect::<Vec<Nt>>()
    );
    println!(
        "{:?}",
        (0..4).map(|idx| ht2.get_base(&idx)).collect::<Vec<Nt>>()
    );
    println!(
        "{:?}",
        (0..4).map(|idx| ht3.get_base(&idx)).collect::<Vec<Nt>>()
    );
    println!(
        "{:?}",
        (0..4).map(|idx| ht4.get_base(&idx)).collect::<Vec<Nt>>()
    );

    println!("---get_sequence---");
    println!("{:?}", wt.get_sequence());
    println!("{:?}", ht.get_sequence());
    println!("{:?}", ht2.get_sequence());
    println!("{:?}", ht3.get_sequence());
    println!("{:?}", ht4.get_sequence());
}
