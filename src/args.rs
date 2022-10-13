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
    /// Number of generations to simulate.
    #[clap(short, long, default_value_t = 200)]
    pub generations: usize,

    /// Path to settings.
    #[clap(long)]
    pub settings: String,

    /// Path to transfer plan.
    #[clap(long)]
    pub transfer_plan: String,

    /// Path to sequence. Format: FASTA.
    #[clap(long)]
    pub sequence: String,

    /// Path to output. Format: FASTA.
    #[clap(long, short)]
    pub output: String,

    /// Path to log file
    #[clap(long, default_value = "virolution.log")]
    pub log_file: String,

    /// Path to optional tree file. Format: Extended Newick.
    #[clap(long)]
    pub trees: Option<String>,

    /// Number of compartments
    #[clap(long, default_value = "3")]
    pub n_compartments: usize,
}
