use crate::index_builder::core::build_index;
use anyhow::Result;
use clap::Parser;
use memmap2::Mmap;
use std::{
    fs::File,
    io::{BufWriter, IoSlice, Write, stdout},
    path::{Path, PathBuf},
    str,
};

pub fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let filename = path.file_name().unwrap_or_default().to_string_lossy();
    parent.join(format!("{filename}{suffix}"))
}

#[derive(Debug, Clone, Parser)]
pub struct CommonArgs {
    /// Input GFF file path
    #[arg(short, long)]
    pub input: PathBuf,

    /// Output file (stdout if not provided)
    #[arg(short, long)]
    pub output: Option<PathBuf>,

    /// Comma-separated feature types to retain (e.g. exon,gene)
    #[arg(short = 'T', long = "types")]
    pub types: Option<String>,

    /// Enable verbose output
    #[arg(short = 'v', long, default_value_t = false)]
    pub verbose: bool,

    /// Number of threads for parallel processing
    #[arg(short = 't', long, default_value_t = 4)]
    pub threads: usize,
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

/// Optimized write_gff_output: memory-maps input, merges adjacent/overlapping blocks,
/// and uses buffered output to minimize syscalls.
pub fn write_gff_output(
    gff_path: &Path,
    blocks: &[(u64, u64)],
    output_path: &Option<std::path::PathBuf>,
    _allowed_types: Option<&str>,
    verbose: bool,
) -> Result<()> {
    // Open and mmap the GFF file
    let file = File::open(gff_path)?;
    let mmap = unsafe { Mmap::map(&file)? };

    // Sort and merge blocks to reduce number of slices
    let mut sorted = blocks.to_vec();
    sorted.sort_unstable_by_key(|&(start, _)| start);
    let mut merged = Vec::with_capacity(sorted.len());
    let mut iter = sorted.into_iter();
    if let Some((mut cur_start, mut cur_end)) = iter.next() {
        for (start, end) in iter {
            if start <= cur_end {
                // overlapping or adjacent
                cur_end = cur_end.max(end);
            } else {
                merged.push((cur_start, cur_end));
                cur_start = start;
                cur_end = end;
            }
        }
        merged.push((cur_start, cur_end));
    }

    // Prepare buffered writer
    let out_file: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(stdout())),
    };
    let mut writer = out_file;

    // Collect IoSlice references
    let mut slices: Vec<IoSlice> = Vec::with_capacity(merged.len());
    for &(so, eo) in &merged {
        let start = so as usize;
        let len = (eo - so) as usize;
        slices.push(IoSlice::new(&mmap[start..start + len]));
    }

    // Write all in one vectored syscall (or minimal syscalls)
    let mut written = 0;
    while written < slices.len() {
        let nw = writer.write_vectored(&slices[written..])?;
        // Determine how many IoSlices were fully written
        let mut consumed = 0;
        let mut remaining = nw;
        for slice in &slices[written..] {
            if remaining >= slice.len() {
                remaining -= slice.len();
                consumed += 1;
            } else {
                break;
            }
        }
        written += consumed;
        if consumed == 0 {
            // fallback to basic write to make progress
            let fallback = &slices[written][..];
            let _ = writer.write(fallback)?;
            written += 1;
        }
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

// Checks if all required index files exist
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
