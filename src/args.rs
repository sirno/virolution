use clap::Parser;

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about,
    long_about = None,
    name = "virolution",
    trailing_var_arg = true,
)]
pub struct Args {
    /// Name of the simulation
    #[clap(short, long, default_value = "simulation")]
    pub name: String,

    /// Number of generations to simulate
    #[clap(short, long, default_value_t = 200)]
    pub generations: usize,

    /// Path to settings. Format: YAML
    #[clap(long, default_value = "settings.yaml")]
    pub settings: String,

    /// Path to sequence. Format: FASTA
    #[clap(long, default_value = "reference.fasta")]
    pub sequence: String,

    /// Path to output directory
    #[clap(long)]
    pub outdir: Option<String>,

    /// Sampling format
    #[clap(long, default_value = "fasta")]
    pub sampling_format: String,

    /// Path to fitness table
    #[clap(long, default_value = "fitness_table_{}.npy")]
    pub fitness_table: String,

    /// Path to log file
    #[clap(long, default_value = "virolution.log")]
    pub logfile: String,

    /// Path to optional tree file. Format: Extended Newick
    #[clap(long)]
    pub trees: Option<String>,

    /// Path to optional ancestry tree file which contains the ancestry of each
    /// individual at the last iteration. Format: Extended Newick
    #[clap(long)]
    pub ancestry: Option<String>,

    /// Number of compartments
    #[clap(long, default_value = "3")]
    pub n_compartments: usize,

    /// Initial population size
    #[clap(long)]
    pub initial_population_size: Option<usize>,

    /// Initial population source file
    #[clap(long)]
    pub initial_population_file: Option<String>,

    /// Disable progress bar
    #[clap(long)]
    pub disable_progress_bar: bool,

    /// Number of threads.
    /// If not set, the number of threads is set to the number of logical cores
    #[clap(long)]
    pub threads: Option<usize>,

    /// Verbosity
    #[clap(long, short, action = clap::ArgAction::Count, default_value = "0")]
    pub verbose: u8,

    /// Test mode for use during development
    #[cfg(debug_assertions)]
    #[clap(long)]
    pub test: bool,
}
