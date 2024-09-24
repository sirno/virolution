#![feature(let_chains)]
#![feature(get_mut_unchecked)]
#![feature(extract_if)]
#![feature(test)]

extern crate test;

use clap::Parser;
use virolution::args::Args;
use virolution::runner::Runner;

fn main() {
    let args = Args::parse();

    if cfg!(feature = "parallel") {
        println!("Running in parallel mode.");
    } else {
        println!("Running in serial mode.");
    }

    let mut runner = Runner::new(args).unwrap_or_else(|err| {
        eprintln!("Unable to init simulation: {err}.");
        std::process::exit(1);
    });

    runner.start();
}

/// TODO: These tests should run for multiple feature flags
#[cfg(test)]
mod tests {
    use std::process::Command;

    struct Cleanup<'a> {
        dir: &'a str,
    }

    impl<'a> Drop for Cleanup<'a> {
        fn drop(&mut self) {
            std::fs::remove_dir_all(self.dir).unwrap();
        }
    }

    /// Test invocation of program for singlehost setup with 1 generation
    #[test]
    fn test_singlehost() {
        let tmp_dir = ".tmp_singlehost";
        let _cleanup = Cleanup { dir: tmp_dir };
        let output = Command::new("cargo")
            .args([
                "run",
                "--",
                "--disable-progress-bar",
                "--settings",
                "test/singlehost/singlehost.yaml",
                "--sequence",
                "test/singlehost/ref.fasta",
                "--generations",
                "1",
                "--n-compartments",
                "2",
                "--outdir",
                tmp_dir,
            ])
            .output()
            .expect("Failed to execute command");

        assert!(output.status.success());

        // Check if output files are created
        assert!(std::path::Path::new(&format!("{}/fitness_table_0.npy", tmp_dir)).exists());
        assert!(std::path::Path::new(&format!("{}/virolution.log", tmp_dir)).exists());
    }

    /// Test invocation of program for multihost setup with 1 generation
    #[test]
    fn test_multihost() {
        let tmp_dir = ".tmp_multihost";
        let _cleanup = Cleanup { dir: tmp_dir };

        let output = Command::new("cargo")
            .args([
                "run",
                "--",
                "--disable-progress-bar",
                "--settings",
                "test/multihost/multihost.yaml",
                "--sequence",
                "test/multihost/ref.fasta",
                "--generations",
                "1",
                "--n-compartments",
                "2",
                "--outdir",
                tmp_dir,
            ])
            .output()
            .expect("Failed to execute command");

        assert!(output.status.success());

        // Check if output files are created
        assert!(std::path::Path::new(&format!("{}/fitness_table_0.npy", tmp_dir)).exists());
        assert!(std::path::Path::new(&format!("{}/fitness_table_1.npy", tmp_dir)).exists());
        assert!(std::path::Path::new(&format!("{}/virolution.log", tmp_dir)).exists());
    }

    /// Test invocation of program for epistasis setup with 1 generation
    #[test]
    fn test_epistasis() {
        let tmp_dir = ".tmp_epistasis";
        let output = Command::new("cargo")
            .args([
                "run",
                "--",
                "--disable-progress-bar",
                "--settings",
                "test/epistasis/epistasis.yaml",
                "--sequence",
                "test/epistasis/ref.fasta",
                "--generations",
                "1",
                "--n-compartments",
                "2",
                "--outdir",
                tmp_dir,
            ])
            .output()
            .expect("Failed to execute command");

        assert!(output.status.success());

        // Check if output files are created
        assert!(std::path::Path::new(&format!("{}/fitness_table_0.npy", tmp_dir)).exists());
        assert!(std::path::Path::new(&format!("{}/epistasis_table_0.npy", tmp_dir)).exists());
        assert!(std::path::Path::new(&format!("{}/virolution.log", tmp_dir)).exists());
    }
}
