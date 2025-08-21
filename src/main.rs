use anyhow::Result;
use clap::{Parser, Subcommand};
use gffx::index;
//use gffx::sort;
use gffx::extract;
use gffx::intersect;
use gffx::search;

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
    Index(index::IndexArgs),
    //    Sort(sort::SortArgs),
    Intersect(intersect::IntersectArgs),
    Extract(extract::ExtractArgs),
    Search(search::SearchArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();
    match cli.command {
        Commands::Index(args) => gffx::index::run(&args)?,
        //       Commands::Sort(args) => gffx::sort::run(&args)?,
        Commands::Intersect(args) => gffx::intersect::run(&args)?,
        Commands::Extract(args) => gffx::extract::run(&args)?,
        Commands::Search(args) => gffx::search::run(&args)?,
    }
    Ok(())
}
