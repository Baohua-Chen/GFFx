use crate::{CommonArgs, load_fts, load_gof, load_prt, write_gff_output};
use anyhow::{Context, Result, bail};
use clap::Parser;
use rayon::prelude::*;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
    sync::{Arc, Mutex},
    time::Instant,
};

/// Extract subtrees from a GFF file by a list of feature names (from --feature-file).
#[derive(Parser, Debug)]
#[command(
    about = "Extract models by feature IDs",
    long_about = "This tool extracts features and their parent models by feature IDs"
)]
#[clap(group(
    clap::ArgGroup::new("feature")
        .required(true)
        .args(&["feature_file", "feature_id"])
))]
pub struct ExtractArgs {
    #[clap(flatten)]
    pub common: CommonArgs,

    #[arg(short = 'f', long, group = "feature")]
    pub feature_id: Option<String>,

    #[arg(short = 'F', long, group = "feature")]
    pub feature_file: Option<PathBuf>,

    #[arg(short = 'd', long)]
    pub descendants_only: bool,
}

/// Recursively walk up parent map to find root feature ID, with loop detection
fn resolve_root(start: u32, prt_map: &HashMap<u32, u32>) -> Result<u32> {
    let mut current = start;
    let mut visited = HashSet::new();
    while visited.insert(current) {
        match prt_map.get(&current) {
            Some(&parent) if parent == current => return Ok(current),
            Some(&parent) => current = parent,
            None => bail!("Parent not found for feature ID {}", current),
        }
    }
    bail!(
        "Circular parent relationship detected for feature ID {}",
        start
    )
}

fn extract_features(
    feature_ids: &HashSet<String>,
    fts_map: &HashMap<&str, u32>,
    prt_map: &HashMap<u32, u32>,
    gof_map: &HashMap<u32, (u64, u64)>,
    descendants_only: bool,
    threads: usize,
) -> Vec<(u64, u64)> {
    let seen_roots: Arc<Mutex<HashSet<u32>>> = Arc::new(Mutex::new(HashSet::new()));
    let process_feature = |fid: &String| -> Option<(u64, u64)> {
        let numeric_id = *fts_map.get(fid.as_str())?;
        let root_id = resolve_root(numeric_id, prt_map).ok()?;

        {
            let mut seen = seen_roots.lock().unwrap();
            if !seen.insert(root_id) {
                return None;
            }
        }

        if descendants_only {
            eprintln!("Descendants only mode will implemented in the next version.");
            return None;
        }

        gof_map.get(&root_id).cloned()
    };
    if threads > 1 && feature_ids.len() > 2 {
        feature_ids.par_iter().filter_map(process_feature).collect()
    } else {
        feature_ids.iter().filter_map(process_feature).collect()
    }
}

pub fn run(args: &ExtractArgs) -> Result<()> {
    let gff_path = &args.common.input;

    // Start overall timer
    let overall_start = Instant::now();
    let verbose = args.common.verbose;
    if verbose {
        eprintln!("[timing] Starting processing of {:?}", gff_path);
    }

    // Build thread pool
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.common.threads)
        .build_global()
        .expect("Failed to set thread pool size");

    // Load features
    let fts = load_fts(gff_path)?;

    // Load parent relations
    let prt = load_prt(gff_path)?;

    // Load GFF offsets

    let gof = load_gof(gff_path)?;

    // Build lookup maps
    let fts_map: HashMap<&str, u32> = fts
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i as u32))
        .collect();
    let prt_map: HashMap<u32, u32> = prt.iter().map(|e| (e.child, e.parent)).collect();
    let gof_map: HashMap<u32, (u64, u64)> = gof
        .iter()
        .map(|e| (e.feature_id, (e.start_offset, e.end_offset)))
        .collect();

    // Read feature IDs
    let feature_ids: HashSet<String> = if let Some(ref file_path) = args.feature_file {
        let file = File::open(file_path)
            .with_context(|| format!("Cannot open feature list: {:?}", file_path))?;
        let reader = BufReader::new(file);

        reader
            .lines() // Iterator<Item = io::Result<String>>
            .try_fold(
                HashSet::new(),
                |mut set, line| -> Result<HashSet<String>, std::io::Error> {
                    let s = line?;
                    let s = s.trim();
                    if !s.is_empty() {
                        set.insert(s.to_owned());
                    }
                    Ok(set)
                },
            )?
    } else if let Some(ref single_id) = args.feature_id {
        [single_id.clone()].into_iter().collect()
    } else {
        bail!("Either --feature-id (-f) or --feature-file (-F) must be specified");
    };

    // Extract features
    let block_offsets = extract_features(
        &feature_ids,
        &fts_map,
        &prt_map,
        &gof_map,
        args.descendants_only,
        args.common.threads,
    );

    // Write output
    write_gff_output(
        gff_path,
        &block_offsets,
        &args.common.output,
        args.common.types.as_deref(),
        args.common.verbose,
    )?;

    if verbose {
        eprintln!("[timing] Total elapsed: {:?}", overall_start.elapsed());
    }

    Ok(())
}
