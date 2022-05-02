mod haplotype;

use haplotype::Wildtype;

fn main() {
    let bytes = vec![0x00, 0x00, 0x00, 0x00];
    let wt = Wildtype::create_wildtype(bytes);
    println!("wt: {:?}", wt);
    let ht = wt.borrow_mut().create_descendant(2, 0x01);
    let ht2 = ht.borrow_mut().create_descendant(1, 0x02);
    let ht3 = ht2.borrow_mut().create_descendant(2, 0x03);
    println!("ht: {:?}", ht);
    println!("ht2: {:?}", ht2);
    println!("ht3: {:?}", ht3);
    for i in 0..4 {
        println!("{}: {}", i, ht3.borrow().get_base(i));
    }
    let seq = ht2.borrow().get_sequence();
    println!("{:?}", seq);
}
