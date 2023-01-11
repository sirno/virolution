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
    /// Number of generations to simulate
    #[clap(short, long, default_value_t = 200)]
    pub generations: usize,

    /// Path to settings. Format: YAML
    #[clap(long)]
    pub settings: String,

    /// Path to transfer plan. Format: CSV
    #[clap(long)]
    pub transfer_plan: String,

    /// Path to sequence. Format: FASTA
    #[clap(long)]
    pub sequence: String,

    /// Path to samples directory
    #[clap(long, default_value = "./samples")]
    pub outdir: String,

    /// Path to fitness table
    #[clap(long, default_value = "fitness_table.npy")]
    pub fitness_table: String,

    /// Simulation name
    #[clap(long, default_value = "simulation")]
    pub simulation_name: String,

    /// Path to log file
    #[clap(long, default_value = "virolution.log")]
    pub log_file: String,

    /// Path to optional tree file. Format: Extended Newick
    #[clap(long)]
    pub trees: Option<String>,

    /// Number of compartments
    #[clap(long, default_value = "3")]
    pub n_compartments: usize,

    /// Initial population size
    #[clap(long, default_value = "10000000")]
    pub initial_population_size: usize,

    /// Disable progress bar
    #[clap(long)]
    pub disable_progress_bar: bool,
}
