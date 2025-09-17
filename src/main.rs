use anyhow::Result;
use clap::{Parser, Subcommand};
use gffx::commands::*;

#[derive(Parser)]
#[command(
    name = "gffx",
    version,
    about = concat!("GFFx: A ultra-fast feature extractor for GFF files\nVersion: ", env!("CARGO_PKG_VERSION")),
    propagate_version = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Index(IndexArgs),
    Intersect(IntersectArgs),
    Extract(ExtractArgs),
    Search(SearchArgs),
    Coverage(CoverageArgs),
    Depth(DepthArgs),
    Sample(SampleArgs)
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index(args) => run_index(&args)?,
        Commands::Intersect(args) => run_intersect(&args)?,
        Commands::Extract(args) => run_extract(&args)?,
        Commands::Search(args) => run_search(&args)?,
        Commands::Coverage(args) => run_coverage(&args)?,
        Commands::Depth(args) => run_depth(&args)?,
        Commands::Sample(args) => run_sample(&args)?,
    }

    Ok(())
}
