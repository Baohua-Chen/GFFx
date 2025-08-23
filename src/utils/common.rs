use crate::index_builder::core::build_index;
use anyhow::{Result, bail};
use clap::Parser;
use memmap2::Mmap;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    fs::File,
    io::{BufWriter, IoSlice, Write, stdout},
    path::{Path, PathBuf},
    str,
};

#[derive(Debug, Clone, Parser)]
pub struct CommonArgs {
    /// Input GFF file path
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Output file (stdout if not provided)
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub output: Option<PathBuf>,

    #[arg(short = 'F', long = "full-model", default_value_t = false)]
    pub full_model: bool,

    /// Comma-separated feature types to retain (e.g. exon,gene)
    #[arg(short = 'T', long = "types", value_name = "TYPES")]
    pub types: Option<String>,

    /// Enable verbose output
    #[arg(
        short = 'v',
        long = "verbose",
        default_value_t = false,
        value_name = "BOOL",
        help = "Enable verbose output"
    )]
    pub verbose: bool,

    /// Number of threads for parallel processing
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 12,
        value_name = "NUM",
        help = "Number of rayon threads (default: num_cpus)"
    )]
    pub threads: usize,
}

impl CommonArgs {
    #[inline]
    pub fn effective_threads(&self) -> usize {
        if self.threads == 0 {
            // Use all available cores when user didn't specify
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        } else {
            self.threads
        }
    }

    pub fn init_rayon(&self) {
        let n = self.effective_threads();

        match rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
        {
            Ok(()) => {
                if self.verbose {
                    eprintln!("[INFO] rayon threads = {}", n);
                }
            }
            Err(e) => {
                if self.verbose {
                    eprintln!("[WARN] rayon global pool already initialized: {e}");
                }
            }
        }
    }
}

/* Append a suffix string to a given file path, returning the new path.

    This is commonly used to generate paths for derived index files
    (e.g., appending `.gof` or `.fts` to a GFF file path).

    # Arguments
    - `path`: The original file path.
    - `suffix`: The string to append to the filename.

    # Returns
    A new `PathBuf` with the suffix added.

    # Example
    ```rust
    # use std::path::Path;
    # use your_crate::append_suffix;
    let p = append_suffix(Path::new("/data/a.gff"), ".gof");
    assert_eq!(p.to_string_lossy(), "/data/a.gff.gof");
*/
pub fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let filename = path.file_name().unwrap_or_default().to_string_lossy();
    parent.join(format!("{filename}{suffix}"))
}

/// Write GFF header lines (starting with '#') to output
/// Returns the byte position after the header
pub fn write_gff_header<W: Write>(writer: &mut W, gff_buf: &[u8]) -> Result<usize> {
    let mut pos = 0;
    while pos < gff_buf.len() && gff_buf[pos] == b'#' {
        if let Some(nl) = gff_buf[pos..].iter().position(|&b| b == b'\n') {
            let end = pos + nl + 1;
            writer.write_all(&gff_buf[pos..end])?;
            pos = end;
        } else {
            break;
        }
    }
    Ok(pos)
}

/// Write selected byte ranges ("blocks") of a GFF file to an output file or stdout.
///
/// Features:
/// - Uses memory-mapped I/O for efficiency.
/// - Merges adjacent or overlapping ranges before writing.
/// - Uses vectored I/O (`write_vectored`) to minimize syscalls.
///
/// # Arguments
/// - `gff_path`: Path to the source GFF file.
/// - `blocks`: A list of `(start, end)` byte ranges to extract.
/// - `output_path`: Output file path. If `None`, writes to stdout.
/// - `_allowed_types`: Reserved for future filtering by feature type (currently unused).
/// - `verbose`: Whether to print diagnostic output.
///
/// # Errors
/// Returns any I/O or mmap errors.
pub fn write_gff_output(
    gff_path: &Path,
    blocks: &[(u64, u64)],
    output_path: &Option<std::path::PathBuf>,
    _allowed_types: Option<&str>,
    verbose: bool,
) -> Result<()> {
    let file = File::open(gff_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let file_len = mmap.len();

    // sort and merge blocks
    let mut sorted = blocks.to_vec();
    sorted.sort_unstable_by_key(|&(s, _)| s);

    let mut merged: Vec<(u64, u64)> = Vec::with_capacity(sorted.len());
    let mut it = sorted.into_iter();
    if let Some((mut cs, mut ce)) = it.next() {
        for (s, e) in it {
            if s <= ce {
                ce = ce.max(e);
            } else {
                if cs < ce {
                    merged.push((cs, ce));
                }
                cs = s;
                ce = e;
            }
        }
        if cs < ce {
            merged.push((cs, ce));
        }
    }

    // build IoSlice list
    let mut slices: Vec<IoSlice<'_>> = Vec::with_capacity(merged.len());
    for &(so, eo) in &merged {
        if so >= eo {
            continue;
        }
        let (start, end) = (so as usize, eo as usize);
        if end > file_len {
            continue;
        }
        slices.push(IoSlice::new(&mmap[start..end]));
    }

    // Write in batches
    let mut writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(stdout())),
    };

    const MAX_IOV: usize = 1024;
    let mut base = 0;
    while base < slices.len() {
        let end = (base + MAX_IOV).min(slices.len());
        let batch = &slices[base..end];

        let nw = writer.write_vectored(batch)?;
        let mut remaining = nw;
        let mut i = 0;

        while i < batch.len() && remaining >= batch[i].len() {
            remaining -= batch[i].len();
            i += 1;
        }

        if i < batch.len() && remaining > 0 {
            let cur = &batch[i];
            writer.write_all(&cur[remaining..])?;
            i += 1;
        }

        for s in &batch[i..] {
            writer.write_all(s)?;
        }

        base = end;
    }

    writer.flush()?;

    if verbose {
        eprintln!(
            "Wrote {} merged GFF block(s) with vectored I/O",
            merged.len()
        );
    }
    Ok(())
}

/// Check if all expected index files for a given GFF exist.
///
/// Expected suffixes: `.gof`, `.fts`, `.prt`, `.sqs`, `.atn`, `.a2f`, `.rit`, `.rix`.
///
/// If any are missing:
/// - If `rebuild = true`, rebuild the index in place and return `Ok(true)`.
/// - Otherwise, return `Ok(false)`.
pub fn check_index_files_exist(
    gff: &PathBuf,
    rebuild: bool,
    attr_key: &str,
    verbose: bool,
) -> Result<bool> {
    let expected_suffixes = [
        ".gof", ".fts", ".prt", ".sqs", ".atn", ".a2f", ".rit", ".rix",
    ];
    let mut missing = Vec::new();

    for ext in &expected_suffixes {
        let path = append_suffix(gff, ext);
        if !path.exists() {
            missing.push(ext.to_string());
        }
    }

    if !missing.is_empty() {
        if rebuild {
            build_index(gff, attr_key, verbose)?;
            return Ok(true);
        }
        return Ok(false);
    }
    Ok(true)
}

/// A grouping of matched feature IDs under their root.
pub struct RootMatched {
    /// Root node numeric ID
    pub root: u32,
    /// Matched numeric IDs that belong under this root
    pub matched: Vec<u32>,
}

/// Resolve the root of a feature ID by following parent links.
///
/// - If the parent points to itself, it's the root.
/// - If no parent is found, treat it as a root.
/// - If a cycle is detected, return an error.
pub fn resolve_root(start: u32, prt_map: &FxHashMap<u32, u32>) -> Result<u32> {
    let mut current = start;
    let mut visited = FxHashSet::default();
    while visited.insert(current) {
        match prt_map.get(&current) {
            Some(&parent) if parent == current => return Ok(current),
            Some(&parent) => current = parent,
            // Treat missing parent as a root as well; adjust if your data guarantees self-parented roots.
            None => return Ok(current),
        }
    }
    bail!(
        "Circular parent relationship detected for feature ID {}",
        start
    )
}

/// Group feature IDs by their resolved roots.
///
/// Converts string feature IDs into numeric IDs, resolves their root
/// via `resolve_root`, and then groups them under each root.
///
/// Uses Rayon parallelism if `threads > 1` and more than 2 features.
///
/// # Returns
/// A `Vec<RootMatched>`, each containing one root and its matched descendants.
pub fn extract_root_matches(
    feature_ids: &FxHashSet<String>,
    fts_map: &FxHashMap<&str, u32>,
    prt_map: &FxHashMap<u32, u32>,
    threads: usize,
) -> Vec<RootMatched> {
    // Fold into HashMap<root, HashSet<numeric_id>>
    let fold_one = |mut acc: FxHashMap<u32, FxHashSet<u32>>, fid: &String| {
        if let Some(&num) = fts_map.get(fid.as_str())
            && let Ok(root) = resolve_root(num, prt_map)
        {
            acc.entry(root).or_default().insert(num);
        }
        acc
    };
    let reduce_maps = |mut a: FxHashMap<u32, FxHashSet<u32>>, b: FxHashMap<u32, FxHashSet<u32>>| {
        for (k, vs) in b {
            a.entry(k).or_default().extend(vs);
        }
        a
    };

    let grouped: FxHashMap<u32, FxHashSet<u32>> = if threads > 1 && feature_ids.len() > 2 {
        feature_ids
            .par_iter()
            .fold(FxHashMap::default, fold_one)
            .reduce(FxHashMap::default, reduce_maps)
    } else {
        feature_ids.iter().fold(FxHashMap::default(), fold_one)
    };

    let mut out = Vec::with_capacity(grouped.len());
    for (root, set) in grouped {
        out.push(RootMatched {
            root,
            matched: set.into_iter().collect(),
        });
    }
    out
}

/// Map root IDs to their byte offsets in the GFF file.
///
/// Looks up each root in a `gof_map` and returns its `(start, end)` offsets.
///
/// # Returns
/// A list of `(root_id, start_offset, end_offset)` tuples.
pub fn roots_to_offsets(
    roots: &[u32],
    gof_map: &FxHashMap<u32, (u64, u64)>,
) -> Vec<(u32, u64, u64)> {
    let mut out = Vec::with_capacity(roots.len());
    for &r in roots {
        if let Some(&(s, e)) = gof_map.get(&r) {
            out.push((r, s, e));
        }
    }
    out
}
