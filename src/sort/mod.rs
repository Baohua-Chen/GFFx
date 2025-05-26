use anyhow::{Result, Context};
use clap::Parser;
use memmap2::MmapOptions;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
};

use crate::{CommonArgs, load_gbi, load_gof, write_gff_output};

/// Command-line arguments for sorting a GFF file by root feature
#[derive(Parser, Debug)]
#[command(
    about = "Sort a GFF file by root feature sequence and position",
    long_about = "This command sorts or reverse-sorts a GFF file's feature tree based on root feature metadata, with optional feature-type filtering."
)]
pub struct SortArgs {
    /// Add common arguments
    #[clap(flatten)]
    pub common: CommonArgs,

    /// Reverse sort order (descending)
    #[arg(short = 'r', long = "reverse", default_value_t = false)]
    pub reverse: bool,
}

/// Holds file offset metadata for a root feature
struct RootInfo {
    seqid_id:  u32,
    start:     u32,
    gof_start: u64,
    gof_end:   u64,
}

/// Main logic: load indexes, sort, and write output
pub fn run(args: &SortArgs) -> Result<()> {
    // Memory-map input GFF
    if args.common.verbose {
        eprintln!("Mapping GFF file: {:?}", args.common.input);
    }
    let file = File::open(&args.common.input)
        .with_context(|| format!("Cannot open GFF file: {:?}", args.common.input))?;
    let mmap = unsafe { MmapOptions::new().map(&file) }
        .with_context(|| format!("Cannot memory-map GFF file: {:?}", args.common.input))?;
    let gff_buf = &mmap[..];

    // Load GOF and GBI indexes
    if args.common.verbose {
        eprintln!("Loading GOF and GBI indexes...");
    }
    let gof_entries = load_gof(&args.common.input).context("Failed to load GOF index")?;
    let gbi_entries = load_gbi(&args.common.input).context("Failed to load GBI index")?;

    // Build feature ID â (start_offset, end_offset) map
    let mut gof_map: HashMap<u32, (u64, u64)> = HashMap::new();
    for entry in gof_entries {
        gof_map.insert(entry.feature_id, (entry.start_offset, entry.end_offset));
    }

    // Deduplicate GBI entries to one root per feature ID
    let mut seen = HashSet::new();
    let mut roots = Vec::new();
    for entry in gbi_entries {
        if seen.insert(entry.feature_id) {
            if let Some(&(so, eo)) = gof_map.get(&entry.feature_id) {
                roots.push(RootInfo {
                    seqid_id:  entry.seqid_id,
                    start:     entry.start,
                    gof_start: so,
                    gof_end:   eo,
                });
            }
        }
    }

    // Sort by seqid_id then start
    if args.common.verbose {
        eprintln!("Sorting {} root features...", roots.len());
    }
    roots.sort_by(|a, b| {
        let ord = a.seqid_id.cmp(&b.seqid_id).then(a.start.cmp(&b.start));
        if args.reverse { ord.reverse() } else { ord }
    });

    // Prepare output writer (file or stdout)
    let blocks: Vec<Vec<u8>> = roots.iter()
        .map(|root| {
            let start = root.gof_start as usize;
            let end = root.gof_end as usize;
            gff_buf[start..end].to_vec()
        })
        .collect();

    write_gff_output(
        gff_buf,
        &blocks,
        &args.common.out,
        args.common.types.as_deref(),
        args.common.verbose,
    )?;

    Ok(())
}

