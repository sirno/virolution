#![feature(get_mut_unchecked)]
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

#[cfg(test)]
mod tests {
    use std::process::Command;

    /// Create a temporary directory for testing and clean it up after the test
    struct TmpDir<'a> {
        dir: &'a str,
    }

    impl<'a> TmpDir<'a> {
        fn new(dir: &'a str) -> Self {
            std::fs::create_dir(dir).unwrap();
            TmpDir { dir }
        }
    }

    impl<'a> Drop for TmpDir<'a> {
        fn drop(&mut self) {
            std::fs::remove_dir_all(self.dir).unwrap();
        }
    }

    #[cfg(feature = "parallel")]
    static RUN_CMD: [&str; 4] = ["run", "--features", "parallel", "--"];

    #[cfg(not(feature = "parallel"))]
    static RUN_CMD: [&str; 2] = ["run", "--"];

    static DEFAULT_ARGS: [&str; 5] = [
        "--disable-progress-bar",
        "--n-compartments",
        "2",
        "--generations",
        "1",
    ];

    /// Test invocation of program for singlehost setup with 1 generation
    #[test]
    fn test_singlehost() {
        let tmp = TmpDir::new(".tmp_singlehost");
        let output = Command::new("cargo")
            .args(RUN_CMD)
            .args(DEFAULT_ARGS)
            .args([
                "--settings",
                "tests/singlehost/singlehost.yaml",
                "--sequence",
                "tests/singlehost/ref.fasta",
                "--outdir",
                tmp.dir,
            ])
            .output()
            .expect("Failed to execute command");

        assert!(output.status.success());

        // Check if output files are created
        assert!(
            std::path::Path::new(&format!("{}/reproductive_fitness_0_table.npy", tmp.dir)).exists()
        );
        assert!(std::path::Path::new(&format!("{}/virolution.log", tmp.dir)).exists());
    }

    /// Test invocation of program for multihost setup with 1 generation
    #[test]
    fn test_multihost() {
        let tmp = TmpDir::new(".tmp_multihost");

        let output = Command::new("cargo")
            .args(RUN_CMD)
            .args(DEFAULT_ARGS)
            .args([
                "--settings",
                "tests/multihost/multihost.yaml",
                "--sequence",
                "tests/multihost/ref.fasta",
                "--outdir",
                tmp.dir,
            ])
            .output()
            .expect("Failed to execute command");

        assert!(output.status.success());

        // Check if output files are created
        assert!(
            std::path::Path::new(&format!("{}/reproductive_fitness_0_table.npy", tmp.dir)).exists()
        );
        assert!(
            std::path::Path::new(&format!("{}/reproductive_fitness_1_table.npy", tmp.dir)).exists()
        );
        assert!(std::path::Path::new(&format!("{}/virolution.log", tmp.dir)).exists());
    }

    /// Test invocation of program for epistasis setup with 1 generation
    #[test]
    fn test_epistasis() {
        let tmp = TmpDir::new(".tmp_epistasis");
        let output = Command::new("cargo")
            .args(RUN_CMD)
            .args(DEFAULT_ARGS)
            .args([
                "--settings",
                "tests/epistasis/epistasis.yaml",
                "--sequence",
                "tests/epistasis/ref.fasta",
                "--outdir",
                tmp.dir,
            ])
            .output()
            .expect("Failed to execute command");

        dbg!(&output);
        assert!(output.status.success());

        // Check if output files are created
        assert!(
            std::path::Path::new(&format!("{}/reproductive_fitness_0_table.npy", tmp.dir)).exists()
        );
        assert!(std::path::Path::new(&format!(
            "{}/reproductive_fitness_0_epistasis_table.npy",
            tmp.dir
        ))
        .exists());
        assert!(std::path::Path::new(&format!("{}/virolution.log", tmp.dir)).exists());
    }
}
