use crate::index_builder::core::build_index;
use anyhow::{Result, Context};
use clap::{Parser, CommandFactory};
use clap::error::ErrorKind;
use memchr::{memchr, memmem};
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

    /// Return the entire gene model for each match (full-model mode); default: only the matched feature (feature-only mode).
    #[arg(short = 'F', long = "full-model", default_value_t = false)]
    pub full_model: bool,

    /// Comma-separated feature types to retain (e.g. exon,gene); only effective in feature-only mode
    #[arg(short = 'T', long = "types", value_name = "TYPES")]
    pub types: Option<String>,

    /// Number of threads for parallel processing
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 12,
        value_name = "NUM",
    )]
    pub threads: usize,

    /// Enable verbose output
    #[arg(
        short = 'v',
        long = "verbose",
        default_value_t = false,
        value_name = "BOOL",
    )]
    pub verbose: bool,
}

impl CommonArgs {
    /// Return the number of effective threads:
    /// - If user sets `--threads 0`, use all available cores
    /// - Otherwise, use the user-specified number
    #[inline]
    pub fn effective_threads(&self) -> usize {
        if self.threads == 0 {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        } else {
            self.threads
        }
    }

    /// Initialize rayon global thread pool
    /// - Uses `effective_threads()` to decide the number of threads
    /// - Prints info/warning if verbose mode is enabled
    pub fn init_rayon(&self) {
        let n = self.effective_threads();
        match rayon::ThreadPoolBuilder::new().num_threads(n).build_global() {
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

    /// Post-parse hook:
    /// - Validate argument combinations
    /// - Print info messages
    /// - Initialize rayon
    pub fn post_parse(&self) -> Result<(), clap::Error> {
        // Custom conflict error: full-model mode cannot be combined with --types
        if self.full_model && self.types.is_some() {
            return Err(
                clap::Error::raw(
                    ErrorKind::ArgumentConflict,
                    "Full-model mode does not support filtering by feature types (-T/--types).",
                )
                .with_cmd(&Self::command())
            );
        }

        // Initialize rayon after validation
        self.init_rayon();

        Ok(())
    }

    /// Combined parse + post-parse helper:
    /// - Parses CLI arguments
    /// - Runs `post_parse()`
    /// - Exits with error if validation fails
    pub fn parse_and_init() -> Self {
        let args = Self::parse();
        if let Err(e) = args.post_parse() {
            e.exit();
        }
        args
    }
}

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
    blocks: &[(u32, u64, u64)],
    output_path: &Option<std::path::PathBuf>,
    verbose: bool,
) -> Result<()> {
    let file = File::open(gff_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let file_len = mmap.len();

    // sort and merge blocks
    let mut sorted: Vec<(u64, u64)> = blocks.iter().map(|&(_, s, e)| (s, e)).collect();
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


pub fn write_gff_output_filtered(
    gff_path: &PathBuf,
    blocks: &[(u32, u64, u64)],
    per_root_matches: &FxHashMap<u32, FxHashSet<String>>,
    atn_attr_name: &str,
    output_path: &Option<PathBuf>,
    types_filter: Option<&str>,
    verbose: bool,
) -> Result<()> {
    // mmap GFF
    let file =
        File::open(gff_path).with_context(|| format!("Cannot open GFF file: {:?}", gff_path))?;
    let mmap =
        unsafe { Mmap::map(&file) }.with_context(|| format!("mmap failed for {:?}", gff_path))?;
    let file_len = mmap.len();

    // parse optional type filter: comma-separated into a HashSet<String>
    let type_allow: Option<FxHashSet<String>> = types_filter.map(|s| {
        s.split(',')
            .map(|t| t.trim().to_string())
            .filter(|t| !t.is_empty())
            .collect()
    });

    let bkey: Vec<u8> = {
        let mut k = atn_attr_name.as_bytes().to_vec();
        k.push(b'=');
        k
    };
    let bkey_finder = memmem::Finder::new(&bkey);
    
    // Process blocks in parallel; each task returns (block_start, matched_bytes)
    let mut parts: Vec<(u64, Vec<u8>)> = blocks
        .par_iter()
        .filter_map(|&(root, start, end)| {
            // root -> set of string IDs to keep
            let keep: &FxHashSet<String> = per_root_matches.get(&root)?;
            if keep.is_empty() {
                return None;
            }

            let s = start as usize;
            let e = end.min(file_len as u64) as usize;
            if s >= e || e > file_len {
                return None;
            }
            let window = &mmap[s..e];

            // Output buffer for this block
            let mut out = Vec::<u8>::with_capacity(1024);
            let mut pos = 0usize;

            // Iterate lines in [s, e)
            let next_line = |from: usize| -> Option<(usize, usize, usize)> {
                if from >= window.len() {
                    return None;
                }
                // '\n' inclusive; rel is the index after '\n' or end-of-window
                let rel = memchr(b'\n', &window[from..])
                    .map(|i| from + i + 1)
                    .unwrap_or(window.len());
                // strip trailing '\n' and optional '\r'
                let mut end_no_nl = rel;
                if end_no_nl > from && window[end_no_nl - 1] == b'\n' {
                    end_no_nl -= 1;
                }
                if end_no_nl > from && window[end_no_nl - 1] == b'\r' {
                    end_no_nl -= 1;
                }
                Some((from, rel, end_no_nl))
            };

            // If a type filter is supplied, check column-3 equals one of allowed types
            let type_ok = |line: &[u8]| -> bool {
                if let Some(allow) = &type_allow {
                    // find first three tabs
                    let i1 = match memchr(b'\t', line) {
                        Some(i) => i,
                        None => return false,
                    };
                    let i2 = match memchr(b'\t', &line[i1 + 1..]) {
                        Some(x) => i1 + 1 + x,
                        None => return false,
                    };
                    let i3 = match memchr(b'\t', &line[i2 + 1..]) {
                        Some(x) => i2 + 1 + x,
                        None => return false,
                    };
                    let ty = &line[i2 + 1..i3];
                    if let Ok(ty_str) = std::str::from_utf8(ty) {
                        allow.contains(ty_str)
                    } else {
                        false
                    }
                } else {
                    true
                }
            };

            // Return true if attributes contain `ID=<value>` and value ∈ keep
            let id_hits_keep = |line_no_crlf: &[u8]| -> bool {
                // move to 9th field (attributes)
                let mut off = 0usize;
                let mut tabs = 0u8;
                while tabs < 8 {
                    match memchr(b'\t', &line_no_crlf[off..]) {
                        Some(i) => {
                            off += i + 1;
                            tabs += 1;
                        }
                        None => return false,
                    }
                }
                let attr = &line_no_crlf[off..];
                if let Some(p) = bkey_finder.find(attr) {
                    let vstart = p + bkey.len();
                    // value ends at ';' or end-of-line
                    let vend = memchr(b';', &attr[vstart..])
                        .map(|i| vstart + i)
                        .unwrap_or(attr.len());
                    let id_slice = &attr[vstart..vend];
                    if let Ok(id_str) = std::str::from_utf8(id_slice) {
                        return keep.contains(id_str);
                    }
                }
                false
            };

            // Scan lines in this block window
            while let Some((ls, le, ln_end)) = next_line(pos) {
                pos = le;
                let line = &window[ls..le];
                if !line.is_empty() && line[0] == b'#' {
                    continue; // skip comments
                }
                let line_no_crlf = &window[ls..ln_end];

                if !type_ok(line_no_crlf) {
                    continue;
                }
                if id_hits_keep(line_no_crlf) {
                    out.extend_from_slice(line);
                }
            }

            if verbose {
                let matched_lines = out.iter().filter(|&&b| b == b'\n').count();
                eprintln!(
                    "[filter] root={} block=[{}..{}] keep_ids={} matched_lines={}",
                    root, start, end, keep.len(), matched_lines
                );
            }

            if out.is_empty() {
                None
            } else {
                Some((start, out))
            }
        })
        .collect();

    // Keep original block order
    parts.sort_unstable_by_key(|(s, _)| *s);

    // Write output (stdout or file)
    let raw: Box<dyn Write> = match output_path {
        Some(p) => Box::new(File::create(p).with_context(|| format!("Cannot create output: {:?}", p))?),
        None => Box::new(std::io::stdout()),
    };
    // Bigger buffer reduces syscalls; tune as needed
    let mut writer = BufWriter::with_capacity(16 * 1024 * 1024, raw);
    for (_, buf) in parts {
        writer.write_all(&buf)?;
    }
    writer.flush()?;
    Ok(())
}


