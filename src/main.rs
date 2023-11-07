#![feature(let_chains)]
#![feature(get_mut_unchecked)]
#![feature(extract_if)]
#![feature(test)]

extern crate test;

use clap::Parser;
use virolution::args::Args;
use virolution::runner::Runner;

fn main() {
    if cfg!(feature = "parallel") {
        println!("Running in parallel mode.");
    } else {
        println!("Running in serial mode.");
    }

    let args = Args::parse();

    let mut runner = Runner::new(args).unwrap_or_else(|err| {
        eprintln!("Unable to init simulation: {err}.");
        std::process::exit(1);
    });

    runner.start();
}
