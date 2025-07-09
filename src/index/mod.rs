use clap::Parser;
use std::path::PathBuf;
use anyhow::Result;
use crate::build_index;

#[derive(Parser, Debug)]
#[command(about = "Build index for GFF file", long_about = "This command builds index files for fast retrieval from a GFF file.")]

pub struct IndexArgs {
    #[arg(short, long)]
    input: PathBuf,

    #[arg(short, long, default_value = "gene_name")]
    pub attribute: String,

    #[arg(short, long, default_value_t = false)]
    verbose: bool,
}

pub fn run(args: &IndexArgs) -> Result<()> {
    if args.verbose {
        println!("Indexing: {}", args.input.display());
    }

    build_index(&args.input, &args.attribute, args.verbose)?;

    if args.verbose {
        println!("Index created successfully.");
    }

    Ok(())
}
