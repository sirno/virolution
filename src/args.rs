use clap::Parser;

#[derive(Parser, Debug)]
#[clap(
    author,
    version,
    about,
    long_about = None,
    name = "value_hints_derive",
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

    /// Path to sequence (fasta file)
    #[clap(long)]
    pub sequence: String,

    /// Path to output (fasta file)
    #[clap(long, short)]
    pub output: String,
}
