use anyhow::{Result, Context};
use std::{
    collections::HashSet, 
    io::{Write, stdout},
    path::{PathBuf, Path},
    fs::File,
    str};
use clap::{Parser};


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
    #[arg(short = 'V', long, default_value_t = false)]
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

/// Writes a GFF block, raw or filtered by the provided comma-separated feature types.
///
/// If `allowed` is `None`, the entire block is written unaltered.
/// Otherwise, only lines whose third column (feature type) is in `allowed` are written,
/// and all comment lines (starting with `#`) are always preserved.
pub fn write_block(
    writer: &mut dyn Write,
    block: &[u8],
    allowed: Option<&str>,
) -> Result<()> {
    if let Some(types_str) = allowed {
        // Build a set of allowed feature types
        let allowed_set: HashSet<String> = types_str
            .split(',')
            .map(|s| s.trim().to_string())
            .collect();
        // Decode block as UTF-8 for line-by-line filtering
        let text = str::from_utf8(block).context("Invalid UTF-8 in GFF block")?;
        for line in text.lines() {
            if line.starts_with('#') {
                // Always write comments
                writeln!(writer, "{}", line)?;
            } else {
                let cols: Vec<&str> = line.split('\t').collect();
                if cols.len() >= 3 && allowed_set.contains(cols[2]) {
                    writeln!(writer, "{}", line)?;
                }
            }
        }
    } else {
        // No filter: write block raw
        writer.write_all(block)?;
    }
    Ok(())
}

pub fn write_gff_output(
    gff_buf: &[u8],
    blocks: &[Vec<u8>],
    output_path: &Option<PathBuf>,
    allowed_types: Option<&str>,
    verbose: bool,
) -> Result<()> {
    let mut writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(File::create(p)?),
        None => Box::new(stdout()),
    };

    let _ = write_gff_header(&mut writer, gff_buf)?;

    for block in blocks {
        write_block(&mut *writer, block, allowed_types)?;
    }

    if verbose {
        eprintln!("Wrote {} GFF block(s)", blocks.len());
    }

    Ok(())
}
