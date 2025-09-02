use crate::{
    CommonArgs,  load_fts, load_gof, load_prt, write_gff_output, write_gff_output_filtered
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
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

    #[arg(short = 'e', long, group = "feature")]
    pub feature_id: Option<String>,

    #[arg(short = 'E', long, group = "feature")]
    pub feature_file: Option<PathBuf>,
}

pub fn run(args: &ExtractArgs) -> Result<()> {
    let gff_path = &args.common.input;

    // Start overall timer
    let overall_start = Instant::now();
    let verbose = args.common.verbose;
    if verbose {
        eprintln!("[DEBUG] Starting processing of {:?}", gff_path);
    }

    // Build thread pool
    args.common.init_rayon();
    if verbose {
        eprintln!(
            "[DEBUG] Thread pool initialized with {} threads",
            args.common.effective_threads()
        );
    }
    // Load features
    let fts = load_fts(gff_path)?;

    // Load parent relations
    let prt = load_prt(gff_path)?;

    // Load GFF offsets
    let gof = load_gof(gff_path)?;

    // Read feature string IDs （feature name）
    let feature_names: FxHashSet<String> = if let Some(ref file_path) = args.feature_file {
        let file = File::open(file_path)
            .with_context(|| format!("Cannot open feature list: {:?}", file_path))?;
        let reader = BufReader::new(file);
        reader.lines().try_fold(
            FxHashSet::default(),
            |mut set, line| -> Result<FxHashSet<String>, std::io::Error> {
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

    // Phase A: group matches by root
    // Phase A: map feature names to numeric fids
    let (fids_set, missing) = fts.map_fnames_to_fids(
        &feature_names,
        args.common.effective_threads()
    );
    if !missing.is_empty() {
        eprintln!("[WARN] {} feature IDs not found: {:?}", missing.len(), missing);
    }

    // Convert set to vec for alignment with roots
    let fid_vec: Vec<u32> = fids_set.iter().copied().collect();

    // Use PrtMap fast resolver to map fid -> root (u32::MAX = invalid)
    let threads = args.common.effective_threads();
    let roots_vec: Vec<u32> = prt.map_fids_to_roots(&fid_vec, threads);

    // Collect invalid fids (print once) and exclude invalid roots
    let mut invalid_fids: Vec<u32> = fid_vec.iter()
        .zip(roots_vec.iter())
        .filter_map(|(&fid, &r)| if r == u32::MAX { Some(fid) } else { None })
        .collect();
    invalid_fids.sort_unstable();
    invalid_fids.dedup();
    if !invalid_fids.is_empty() {
        eprintln!(
            "[WARN] {} numeric feature IDs are invalid (out-of-range child or parent), skipped: {:?}",
            invalid_fids.len(), invalid_fids
        );
    }

    // Deduplicate valid roots (exclude u32::MAX)
    let mut roots: Vec<u32> = roots_vec.iter().copied().filter(|&r| r != u32::MAX).collect();
    roots.sort_unstable();
    roots.dedup();

    // Phase B: roots -> block offsets
    let blocks: Vec<(u32, u64, u64)> = gof.roots_to_offsets(&roots, args.common.effective_threads());

    if !args.common.full_model || args.common.types.is_some() {
        // Build per_root_matches: root_id -> set of STRING feature IDs
        let mut per_root_matches: FxHashMap<u32, FxHashSet<String>> = FxHashMap::default();
        per_root_matches.reserve(roots.len());
        
        // roots_vec[i] is the root, fid_vec[i] is the numeric fid
        for (i, &root) in roots_vec.iter().enumerate() {
            if root == u32::MAX {
                continue;
            }
            // Convert numeric fid -> string ID only once here
            if let Some(id_str) = fts.ids.get(fid_vec[i] as usize) {
                per_root_matches.entry(root).or_default().insert(id_str.clone());
            }
        }
        
        // Emit only exactly matched lines within blocks
        write_gff_output_filtered(
            gff_path,
            &blocks,
            &per_root_matches,
            "ID",
            &args.common.output,
            args.common.types.as_deref(),
            verbose,
        )?;
    } else {
        // Full-model mode: emit full blocks without filtering
        write_gff_output(
            gff_path,
            &blocks,
            &args.common.output,
            verbose,
        )?;
    }

    if verbose {
        eprintln!("[timing] Total elapsed: {:?}", overall_start.elapsed());
    }

    Ok(())
}
