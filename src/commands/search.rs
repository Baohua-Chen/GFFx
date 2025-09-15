use anyhow::{Result, bail};
use clap::{ArgGroup, Parser};
use rustc_hash::{FxHashMap, FxHashSet};
use regex::Regex;
use std::{
    fs::File,
    io::{BufReader, BufRead},
    path::PathBuf,
    time::Instant,
};


use crate::{
    CommonArgs, load_gof, load_prt, load_a2f, load_atn,
    write_gff_output, write_gff_output_filtered,
};

#[derive(Parser, Debug)]
#[command(
    about = "Search features by attribute values",
    group = ArgGroup::new("attr_input")
        .required(true)
        .args(["attr_list", "attr"])
)]
pub struct SearchArgs {
    /// Common input/output/thread arguments
    #[clap(flatten)]
    pub common: CommonArgs,

    #[arg(
        short = 'A',
        long,
        help = "Attribute list file (one per line)",
        group = "attr_input"
    )]
    attr_list: Option<PathBuf>,

    #[arg(
        short = 'a',
        long,
        help = "Single attribute value to search",
        group = "attr_input"
    )]
    attr: Option<String>,

    #[arg(
        short = 'r',
        long,
        help = "Enable regex mode for attribute matching")]
    regex: bool,
}

/// When in `per-feature` mode: within each root block, emit only lines whose `ID` exactly matches
/// the user-specified features under that root. Optional `types_filter` is applied to column 3.
pub fn run(args: &SearchArgs) -> Result<()> {
    let verbose = args.common.verbose;
    let gff_path = &args.common.input;

    // Init thread pool
    let overall_start = Instant::now();
    if verbose {
        eprintln!("[DEBUG] Starting processing of {:?}", gff_path);
        eprintln!(
            "[DEBUG] Thread pool initialized with {} threads",
            args.common.effective_threads()
        );
    }

    // Load index artifacts
    let prt = load_prt(gff_path)?;          // parent pointers (fid -> parent fid)
    let gof = load_gof(gff_path)?;          // GOF offsets (fid -> (start,end))
    let a2f = load_a2f(gff_path)?;          // attribute index -> fid
    let (atn_attr_name, atn_values) = load_atn(gff_path)?; // attribute values table (index-aligned)

    // Collect attribute values from file or single arg
    let attr_values: Vec<String> = if let Some(file) = &args.attr_list {
        let reader = BufReader::new(File::open(file)?);
        reader
            .lines()
            .map(|r| r.map(|s| s.trim().to_owned()))
            .filter(|r| r.as_ref().map_or(true, |s| !s.is_empty()))
            .collect::<Result<Vec<_>, _>>()?
    } else if let Some(val) = &args.attr {
        vec![val.clone()]
    } else {
        bail!("Either --attr-list (-A) or --attr (-a) must be provided.");
    };

    // Step 1: build attribute -> AID list
    // In regex mode, match by regex; otherwise exact string match.
    let mut attr_to_aids: FxHashMap<String, Vec<u32>> = FxHashMap::default();
    if args.regex {
        let patterns: Vec<Regex> = attr_values
            .iter()
            .map(String::as_str)
            .map(Regex::new)
            .collect::<std::result::Result<Vec<_>, _>>()?;

        for (i, val) in atn_values.iter().enumerate() {
            if patterns.iter().any(|re| re.is_match(val)) {
                attr_to_aids.entry(val.clone()).or_default().push(i as u32);
            }
        }
    } else {
        let wanted: FxHashSet<&str> = attr_values.iter().map(String::as_str).collect();
        for (i, val) in atn_values.iter().enumerate() {
            if wanted.contains(val.as_str()) {
                attr_to_aids.entry(val.clone()).or_default().push(i as u32);
            }
        }
    }

    // Nothing matched → early exit with a helpful error
    if attr_to_aids.is_empty() {
        bail!("None of the attributes matched.");
    }

    if verbose {
        eprintln!("[DEBUG] Matched attribute -> AIDs:");
        for (attr_val, aids) in &attr_to_aids {
            eprintln!("  {} => {:?}", attr_val, aids);
        }
    }

    // Step 2: map AIDs -> FIDs via a2f (attribute index to feature id)
    // Note: a2f is expected to be indexable by AID (usize).
    // We also deduplicate per attribute to keep vectors lean.
    let mut attr_to_fids: FxHashMap<String, Vec<u32>> = FxHashMap::default();
    
    for (attr_val, aids) in &attr_to_aids {
        let mut fids = a2f.map_aids_to_fids_vec(aids);
        fids.sort_unstable();
        fids.dedup();
    
        if !fids.is_empty() {
            attr_to_fids.insert(attr_val.clone(), fids);
        }
    }

    if attr_to_fids.is_empty() {
        bail!("No feature IDs (FIDs) resolved from matched attributes.");
    }

    if verbose {
        eprintln!("[DEBUG] Attribute -> FIDs after a2f mapping:");
        for (attr_val, fids) in &attr_to_fids {
            eprintln!("  {} => {:?} ", attr_val, fids);
        }
    }

    // Step 3: map FIDs -> root FIDs using PrtMap::map_fids_to_roots (fast)
    let mut fid_vec: Vec<u32> = attr_to_fids
        .values()
        .flat_map(|v| v.iter().copied())
        .collect();
    fid_vec.sort_unstable();
    fid_vec.dedup();

    if verbose {
        eprintln!("[DEBUG] Total unique FIDs: {}", fid_vec.len());
    }

    let threads = args.common.effective_threads();
    let root: Vec<u32> = prt.map_fids_to_roots(&fid_vec, threads);

    // Collect invalid fids (mapped to u32::MAX), and build a unique root list
    let mut invalid_fids: Vec<u32> = Vec::new();
    let mut roots_effective: Vec<u32> = Vec::with_capacity(root.len());
    for (&fid, r) in fid_vec.iter().zip(root.iter()) {
        if *r == u32::MAX {
            invalid_fids.push(fid);
        } else {
            roots_effective.push(*r);
        }
    }

    if !invalid_fids.is_empty() {
        invalid_fids.sort_unstable();
        invalid_fids.dedup();
        eprintln!(
            "[WARN] {} FIDs have invalid parent chains (or out-of-range): {:?}",
            invalid_fids.len(),
            invalid_fids
        );
    }
    roots_effective.sort_unstable();
    roots_effective.dedup();

    if roots_effective.is_empty() {
        bail!("No valid root features resolved from matched attributes.");
    }
    if verbose {
        eprintln!("[DEBUG] Total unique roots: {}", roots_effective.len());
    }

    let blocks: Vec<(u32, u64, u64)> = gof.roots_to_offsets(&roots_effective, args.common.effective_threads());
    
    if !args.common.entire_group|| args.common.types.is_some() {
        let allowed_roots: FxHashSet<u32> = roots_effective.iter().copied().collect();
        
        let mut fid_to_root: FxHashMap<u32, u32> = FxHashMap::default();
        fid_to_root.reserve(fid_vec.len());
        for (fid, r) in fid_vec.iter().copied().zip(root.iter().copied()) {
            if r != u32::MAX && allowed_roots.contains(&r) {
                fid_to_root.insert(fid, r);
            }
        }
        
        let mut per_root_matches: FxHashMap<u32, FxHashSet<String>> = FxHashMap::default();
        per_root_matches.reserve(roots_effective.len());
        
        for (attr_val, fids) in &attr_to_fids {
            for &fid in fids {
                if let Some(&r) = fid_to_root.get(&fid) {
                    per_root_matches.entry(r).or_default().insert(attr_val.clone());
                }
            }
        }

        write_gff_output_filtered(
            gff_path,
            &blocks,
            &per_root_matches,
            &atn_attr_name,
            &args.common.output,
            args.common.types.as_deref(),
            verbose,
        )?;
    } else {
        let mut fid_to_root: FxHashMap<u32, u32> = FxHashMap::default();
        for (fid, r) in fid_vec.iter().copied().zip(root.clone().into_iter()) {
            if r != u32::MAX {
                fid_to_root.insert(fid, r);
            }
        }
    
        // Step 4: map roots -> (start, end) offsets from GOF
        // Use a cached index to avoid rebuilding a HashMap on every call.
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
