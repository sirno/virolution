// use crate::args::Args;
use clap::CommandFactory;
use clap_complete::{generate_to, shells::Fish};
use std::env;
use std::io::Error;
use std::path::Path;

include!("src/args.rs");

fn main() -> Result<(), Error> {
    let mut cmd = Args::command();
    let home_dir = match env::var_os("HOME") {
        None => return Ok(()),
        Some(d) => d,
    };

    // Fish completions, if fish is installed for user.
    let fish_completions_dir = Path::new(&home_dir).join(".config/fish/completions");
    if fish_completions_dir.exists() {
        let fish_command_path = generate_to(Fish, &mut cmd, "virolution", fish_completions_dir)?;
        println!("cargo:warning=completion file is generated for fish: {fish_command_path:?}",);
    }

    Ok(())
}
