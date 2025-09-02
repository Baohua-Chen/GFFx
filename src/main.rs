use anyhow::Result;
use clap::{Parser, Subcommand};

use gffx::index;
use gffx::intersect;
use gffx::extract;
use gffx::search;

#[derive(Parser)]
#[command(
    name = "gffx",
    version,
    about = concat!("GFFx: An ultra-fast feature extractor for GFF files\nVersion: ", env!("CARGO_PKG_VERSION")),
    propagate_version = true
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Build index for fast queries
    Index(index::IndexArgs),

    /// Interval set intersection on GFF features
    Intersect(intersect::IntersectArgs),

    /// Extract models by feature IDs
    Extract(extract::ExtractArgs),

    /// Search features by patterns
    Search(search::SearchArgs),
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.command {
        Commands::Index(args) => {
            index::run(&args)?
        }
        Commands::Intersect(args) => {
            args.common.post_parse().unwrap_or_else(|e| e.exit());
            intersect::run(&args)?
        }
        Commands::Extract(args) => {
            args.common.post_parse().unwrap_or_else(|e| e.exit());
            extract::run(&args)?
        }
        Commands::Search(args) => {
            args.common.post_parse().unwrap_or_else(|e| e.exit());
            search::run(&args)?
        }
    }

    Ok(())
}
