use anyhow::{Result, bail, Context};
use rayon::prelude::*;
use memmap2::Mmap;
use rust_htslib::bam::{self, Read};
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    str,
    fs::File,
    path::{Path, PathBuf},
    io::{BufWriter, Write},
};
use clap::Parser;
use crate::{
    Interval, TreeIndexData, load_gof, GofMap,
};
use std::time::{Instant, Duration};
use rust_htslib::bam::ext::BamRecordExtensions;

const MISSING: u64 = u64::MAX; // Sentinel for missing entries

// Number of IoSlices per batch writer
const WRITE_BUF_SIZE: usize = 32 * 1024 * 1024;

/// Compute bin index for coordinate `x` given shift k (unused in the new pipeline).
#[inline]
fn _bin_of(x: u32, shift: u32) -> u32 {
    let _ = shift;
    x >> 0
}

/// Arguments
#[derive(Parser, Debug)]
#[command(
    about = "Compute coverage breadth across genomic feature.",
    long_about = "This tool computes sequencing coverage breadth and fraction from high-throughput sequencing (HTS) alignment files (SAM/BAM/CRAM) or user-specified genomic intervals (BED)."
)]
pub struct CoverageArgs {
    /// GFF file path (indexed via GOF)
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Source: BAM/SAM/CRAM or BED
    #[arg(short = 's', long)]
    pub source: PathBuf,
    
    /// Output file (required)
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Number of threads
    #[arg(short = 't', long = "threads", default_value_t = 12, value_name = "NUM")]
    pub threads: usize,

    /// Verbose logs
    #[arg(short = 'v', long = "verbose", default_value_t = false, value_name = "BOOL")]
    pub verbose: bool,
}

/// Fast u32 parse
#[inline(always)]
fn parse_u32_fast(s: &str) -> Option<u32> {
    if s.is_empty() { return None; }
    let mut n: u32 = 0;
    for b in s.as_bytes() {
        let d = b.wrapping_sub(b'0');
        if d > 9 { return None; }
        n = n.checked_mul(10)?.checked_add(d as u32)?;
    }
    Some(n)
}

/// Extract `ID=` from GFF attributes quickly
#[inline(always)]
fn fast_id(attrs: &str) -> Option<&str> {
    let bytes = attrs.as_bytes();
    let mut i = 0;
    while i + 2 < bytes.len() {
        if bytes[i] == b'I' && bytes[i+1] == b'D' && bytes[i+2] == b'=' {
            let mut j = i + 3;
            while j < bytes.len() && bytes[j] != b';' && bytes[j] != b' ' && bytes[j] != b'\t' {
                j += 1;
            }
            return std::str::from_utf8(&bytes[i+3..j]).ok();
        }
        i += 1;
    }
    None
}

/// Disjoint union of intervals assumed to be half-open [s, e)
/// Input may be unsorted and overlapping; output is sorted, non-overlapping.
fn merge_intervals(mut ivs: Vec<(u32,u32)>) -> Vec<(u32,u32)> {
    if ivs.is_empty() { return ivs; }
    ivs.sort_unstable_by_key(|x| x.0);
    let mut out: Vec<(u32,u32)> = Vec::with_capacity(ivs.len());
    let (mut cs, mut ce) = ivs[0];
    for (s,e) in ivs.into_iter().skip(1) {
        if s <= ce {
            if e > ce { ce = e; }
        } else {
            out.push((cs, ce));
            cs = s; ce = e;
        }
    }
    out.push((cs, ce));
    out
}

/// Union length of intervals (half-open). Input may be unsorted.
fn union_len(mut ivs: Vec<(u32,u32)>) -> usize {
    if ivs.is_empty() { return 0; }
    ivs.sort_unstable_by_key(|x| x.0);
    let mut total: usize = 0;
    let (mut cs, mut ce) = ivs[0];
    for (s, e) in ivs.into_iter().skip(1) {
        if s <= ce { ce = ce.max(e); }
        else { total += (ce - cs) as usize; cs = s; ce = e; }
    }
    total += (ce - cs) as usize;
    total
}

/// Collect coverage intervals per root (from BAM/SAM/CRAM).
/// We DO NOT read GFF slices here; only group regions by root_fid.
fn collect_by_root_from_bam(
    bam_path: &Path,
    index_data: &TreeIndexData,
    verbose: bool,
    threads: usize,
) -> Result<FxHashMap<u32, Vec<(u32,u32)>>> {
    let t_open = Instant::now();
    let mut reader = bam::Reader::from_path(bam_path)?;
    reader.set_threads(std::cmp::max(2, threads))?;
    let t_open_elapsed = t_open.elapsed();

    let header = reader.header().to_owned();

    // Build tid -> chr_id mapping
    let t_map_build = Instant::now();
    let mut tid2num: Vec<Option<u32>> = Vec::with_capacity(header.target_count() as usize);
    for tid in 0..header.target_count() {
        let chrom_bytes = header.tid2name(tid).to_owned();
        let chrom = std::str::from_utf8(&chrom_bytes)?;
        let chr_id = index_data.seqid_to_num.get(chrom).copied();
        tid2num.push(chr_id);
    }
    let t_map_build_elapsed = t_map_build.elapsed();

    // Storage: root_fid -> list of raw intervals (to be merged later)
    let mut by_root: FxHashMap<u32, Vec<(u32,u32)>> = FxHashMap::default();

    let mut t_parse = Duration::ZERO;
    let mut t_tidmap = Duration::ZERO;
    let t_tree = Duration::ZERO;

    let mut hits: Vec<&Interval<u32>> = Vec::new();
    for r in reader.records() {
        let t0 = Instant::now();
        let rec = r?;
        if rec.is_unmapped() {
            t_parse += t0.elapsed();
            continue;
        }
        t_parse += t0.elapsed();

        let t1 = Instant::now();
        let tid = rec.tid();
        if tid < 0 {
            t_tidmap += t1.elapsed();
            continue;
        }
        if let Some(chr_id) = tid2num[tid as usize] {
            let start_i64 = rec.pos();
            let end_i64 = rec.reference_end();
            if start_i64 >= 0 && end_i64 > start_i64 {
                let start = (start_i64 as i128).clamp(0, u32::MAX as i128) as u32;
                let end   = (end_i64 as i128).clamp(0, u32::MAX as i128) as u32;

                // Query candidate roots for this region
                if let Some(tree) = index_data.chr_entries.get(&chr_id) {
                    hits.clear();
                    tree.query_interval(start, end, &mut hits);
                    // De-duplicate roots within a single region
                    let mut seen_in_region: FxHashSet<u32> = FxHashSet::default();
                    for h in &hits {
                        if seen_in_region.insert(h.root_fid) {
                            by_root.entry(h.root_fid).or_default().push((start, end));
                        }
                    }
                }
            }
        }
        t_tidmap += t1.elapsed();
    }

    if verbose {
        eprintln!("[TIMER] (1) Opening BAM:       {:.2?}", t_open_elapsed);
        eprintln!("[TIMER] (1b) Build tid2num:    {:.2?}", t_map_build_elapsed);
        eprintln!("[TIMER] (2) Parse+map records: {:.2?}", t_parse + t_tidmap);
        eprintln!("[TIMER] (3) Tree queries:      {:.2?}", t_tree);
        eprintln!("[INFO] Collected {} roots with coverage", by_root.len());
    }

    Ok(by_root)
}

/// Collect coverage intervals per root (from BED).
fn collect_by_root_from_bed(
    bed_path: &Path,
    index_data: &TreeIndexData,
    verbose: bool,
) -> Result<FxHashMap<u32, Vec<(u32,u32)>>> {
    // mmap the entire BED file
    let file = File::open(bed_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = &mmap[..];

    if verbose {
        eprintln!("[INFO] mmap BED file: {} bytes", data.len());
    }

    // collect line offsets
    let mut line_offsets = Vec::with_capacity(1_000_000);
    line_offsets.push(0usize);
    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            line_offsets.push(i + 1);
        }
    }
    if *line_offsets.last().unwrap_or(&0) != data.len() {
        line_offsets.push(data.len());
    }

    let mut by_root: FxHashMap<u32, Vec<(u32,u32)>> = FxHashMap::default();
    let mut hits: Vec<&Interval<u32>> = Vec::new();

    for w in line_offsets.windows(2) {
        let start = w[0];
        let end = w[1];
        if start >= end { continue; }
        let line = &data[start..end];
        if line.is_empty() || line[0] == b'#' { continue; }

        let fields: Vec<&[u8]> = line
            .split(|&b| b == b'\t' || b == b' ')
            .filter(|f| !f.is_empty())
            .collect();
        if fields.len() < 3 { continue; }

        let chrom = match std::str::from_utf8(fields[0]) { Ok(s) => s, Err(_) => continue };
        let s = match std::str::from_utf8(fields[1]).ok().and_then(|x| x.parse::<u32>().ok()) { Some(v) => v, None => continue };
        let e = match std::str::from_utf8(fields[2]).ok().and_then(|x| x.trim_end().parse::<u32>().ok()) { Some(v) => v, None => continue };
        if s >= e { continue; }

        let Some(&chr_num) = index_data.seqid_to_num.get(chrom) else { continue };

        if let Some(tree) = index_data.chr_entries.get(&chr_num) {
            hits.clear();
            tree.query_interval(s, e, &mut hits);
            let mut seen_in_region: FxHashSet<u32> = FxHashSet::default();
            for h in &hits {
                if seen_in_region.insert(h.root_fid) {
                    by_root.entry(h.root_fid).or_default().push((s, e));
                }
            }
        }
    }

    if verbose {
        eprintln!("[INFO] Collected {} roots with coverage (BED)", by_root.len());
    }

    Ok(by_root)
}

/// Compute breadth for all features within a root using pre-merged disjoint coverage.
fn compute_breadth_for_root(
    gff_slice: &[u8],
    cov_merged: &[(u32,u32)], // sorted, non-overlapping coverage intervals
) -> FxHashMap<String, (String, u32, u32, usize)> {
    #[derive(Clone)]
    struct FeatLine {
        id_idx: u32,
        start0: u32, // 0-based inclusive
        end0:   u32, // 0-based exclusive
    }

    // Map ID string -> index
    let mut id_to_idx: FxHashMap<String, u32> = FxHashMap::default();
    let mut id_strings: Vec<String> = Vec::new();
    let mut id_chrom: Vec<String> = Vec::new();

    // Collect feature lines (within this root)
    let mut lines: Vec<FeatLine> = Vec::new();

    if let Ok(text) = str::from_utf8(gff_slice) {
        for line in text.split_terminator('\n') {
            if line.is_empty() || line.as_bytes()[0] == b'#' { continue; }
            let mut cols = line.splitn(9, '\t');
            let (Some(seqid), _, _, Some(start_s), Some(end_s), _, _, _, Some(attrs)) = (
                cols.next(), cols.next(), cols.next(),
                cols.next(), cols.next(), cols.next(),
                cols.next(), cols.next(), cols.next()
            ) else { continue; };

            let (Some(s1), Some(e1)) = (parse_u32_fast(start_s), parse_u32_fast(end_s)) else { continue; };
            if e1 == 0 { continue; }
            let (s1, e1) = if s1 > e1 { (e1, s1) } else { (s1, e1) };
            // GFF is 1-based, inclusive end; convert to half-open 0-based
            let fstart0 = s1.saturating_sub(1);
            let fend0   = e1;

            if let Some(id) = fast_id(attrs) {
                let idx = *id_to_idx.entry(id.to_owned()).or_insert_with(|| {
                    let k = id_strings.len() as u32;
                    id_strings.push(id.to_owned());
                    id_chrom.push(seqid.to_owned());
                    k
                });
                lines.push(FeatLine { id_idx: idx, start0: fstart0, end0: fend0 });
            }
        }
    }

    if lines.is_empty() || cov_merged.is_empty() {
        return FxHashMap::default();
    }

    // Sort features by start for two-pointer sweep against cov_merged
    lines.sort_unstable_by_key(|x| x.start0);

    let mut min_s: Vec<u32> = vec![u32::MAX; id_strings.len()];
    let mut max_e: Vec<u32> = vec![0; id_strings.len()];
    // For each ID, collect overlaps with coverage (we'll union at the end per ID)
    let mut id_overlap: Vec<Vec<(u32,u32)>> = vec![Vec::new(); id_strings.len()];

    // Two-pointer scan: iterate features in order, and advance cov pointer monotonically
    let mut j = 0usize;
    for fl in &lines {
        // Advance coverage pointer until cov[j].end <= feature.start
        while j < cov_merged.len() && cov_merged[j].1 <= fl.start0 {
            j += 1;
        }
        // Record feature span extents for this ID
        if fl.start0 < min_s[fl.id_idx as usize] { min_s[fl.id_idx as usize] = fl.start0; }
        if fl.end0   > max_e[fl.id_idx as usize] { max_e[fl.id_idx as usize] = fl.end0; }

        // Walk through all coverage intervals that might overlap this feature
        let mut k = j;
        while k < cov_merged.len() && cov_merged[k].0 < fl.end0 {
            let s = fl.start0.max(cov_merged[k].0);
            let e = fl.end0.min(cov_merged[k].1);
            if e > s {
                id_overlap[fl.id_idx as usize].push((s, e));
            }
            if cov_merged[k].1 <= fl.end0 {
                k += 1;
            } else {
                break;
            }
        }
    }

    // Finalize per-ID breadth (union of all overlap pieces), and produce outputs
    let mut out: FxHashMap<String, (String, u32, u32, usize)> = FxHashMap::default();
    for i in 0..id_strings.len() {
        let length = if max_e[i] > min_s[i] { (max_e[i] - min_s[i]) as usize } else { 0 };
        let breadth = union_len(std::mem::take(&mut id_overlap[i]));
        if length > 0 || breadth > 0 {
            out.insert(
                id_strings[i].clone(),
                (id_chrom[i].clone(), min_s[i], max_e[i], breadth),
            );
        }
    }
    out
}

/// After collecting raw intervals per root:
/// 1) Merge (union) them into disjoint intervals;
/// 2) Parse GFF slice for that root;
/// 3) Compute breadth/fraction for each feature under this root.
fn finalize_compute_breadth(
    by_root_raw: FxHashMap<u32, Vec<(u32,u32)>>,
    gof: &GofMap,
    gff_mmap: &Mmap,
    threads: usize,
    verbose: bool,
) -> Result<FxHashMap<String, (String, u32, u32, usize)>> {
    let gff_bytes: &[u8] = &gff_mmap[..];
    let idx = gof.index_cached();

    if by_root_raw.is_empty() {
        return Ok(FxHashMap::default());
    }

    let roots_iter = by_root_raw.into_iter();

    // Parallel per-root processing if threads > 1
    let partials: Vec<FxHashMap<String, (String, u32, u32, usize)>> = if threads > 1 {
        roots_iter.par_bridge().map(|(root, ivs)| {
            // Merge coverage intervals for this root
            let cov = merge_intervals(ivs);
            // Locate GFF slice for this root
            match idx.get(&root) {
                Some(&(s_off, e_off)) if s_off != MISSING && e_off != MISSING && e_off > s_off => {
                    let su = usize::try_from(s_off).unwrap();
                    let eu = usize::try_from(e_off).unwrap();
                    compute_breadth_for_root(&gff_bytes[su..eu], &cov)
                }
                _ => FxHashMap::default(),
            }
        }).collect()
    } else {
        let mut v = Vec::new();
        for (root, ivs) in roots_iter {
            let cov = merge_intervals(ivs);
            match idx.get(&root) {
                Some(&(s_off, e_off)) if s_off != MISSING && e_off != MISSING && e_off > s_off => {
                    let su = usize::try_from(s_off).unwrap();
                    let eu = usize::try_from(e_off).unwrap();
                    v.push(compute_breadth_for_root(&gff_bytes[su..eu], &cov));
                }
                _ => v.push(FxHashMap::default()),
            }
        }
        v
    };

    // Merge per-root maps into global results
    let mut global: FxHashMap<String, (String, u32, u32, usize)> = FxHashMap::default();
    for m in partials {
        for (id, (chrom, s, e, b)) in m {
            global.entry(id).and_modify(|(c0, s0, e0, breadth)| {
                if s < *s0 { *s0 = s; }
                if e > *e0 { *e0 = e; }
                // Note: If the same ID appears under multiple roots (rare), breadth is summed.
                // In well-formed GFF partitioning, one ID should belong to a single root.
                *breadth += b;
                let _ = c0;
            }).or_insert((chrom, s, e, b));
        }
    }

    if verbose {
        eprintln!("[INFO] Aggregated {} feature IDs", global.len());
    }

    Ok(global)
}

/// Write "id\tchr\tstart\tend\tbreadth\tfraction" per line.
pub fn write_breadth_results<W: Write>(
    id_map: FxHashMap<String, (String, u32, u32, usize)>,
    mut out: W,
    verbose: bool,
) -> Result<()> {
    use std::fmt::Write as FmtWrite;
    let mut buf = String::with_capacity(WRITE_BUF_SIZE);
    let mut written = 0usize;

    // Header: 6 columns
    writeln!(buf, "id\tchr\tstart\tend\tbreadth\tfraction")?;
    
    for (id, (chr, start, end, breadth)) in id_map {
        let length = end.saturating_sub(start) as usize;
        let fraction = if length > 0 {
            breadth as f64 / length as f64
        } else {
            0.0
        };
        writeln!(buf, "{id}\t{chr}\t{start}\t{end}\t{breadth}\t{:.6}", fraction)?;
        written += 1;

        if buf.len() >= WRITE_BUF_SIZE {
            out.write_all(buf.as_bytes())?;
            buf.clear();
        }
    }

    if !buf.is_empty() {
        out.write_all(buf.as_bytes())?;
    }
    out.flush()?;

    if verbose {
        eprintln!("[INFO] Wrote {written} feature coverage rows.");
    }
    Ok(())
}

/// Main
pub fn run(args: &CoverageArgs) -> Result<()> {
    let verbose = args.verbose;
    let threads = if args.threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        args.threads
    };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global();
    let gff_path = &args.input;
    
    // Step 1: load GOF index + mmap GFF
    let t0 = Instant::now();
    let gof = load_gof(&gff_path)?;
    let file = File::open(gff_path).with_context(|| format!("Cannot open GFF file: {:?}", gff_path))?;
    let gff_mmap = unsafe { Mmap::map(&file) }.with_context(|| format!("GFF mmap failed for {:?}", gff_path))?;
    let t_load_fts = t0.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 1: Load GOF & mmap GFF: {:.2?}", t_load_fts);
    }

    // Step 2: load interval tree index
    let t1 = Instant::now();
    let index_data = TreeIndexData::load_tree_index(&gff_path)?;
    let t_build_index = t1.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 2: Load tree index: {:.2?}", t_build_index);
    }

    // Step 3: collect coverage intervals per root
    let t2 = Instant::now();
    let source_path = &args.source;

    let ext = source_path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_lowercase());

    let by_root = match ext.as_deref() {
        Some("bam") | Some("sam") | Some("cram") => {
            collect_by_root_from_bam(source_path.as_path(), &index_data, verbose, threads)?
        }
        Some("bed") => {
            collect_by_root_from_bed(source_path.as_path(), &index_data, verbose)?
        }
        _ => {
            bail!(
                "Unsupported file type: {:?}. Expected .bam/.sam/.cram or .bed",
                source_path
            );
        }
    };
    let t_collect = t2.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 3: Collect intervals: {:.2?}", t_collect);
    }

    // Step 4: per-root merge & compute breadth over GFF slices
    let t3 = Instant::now();
    let id_map = finalize_compute_breadth(by_root, &gof, &gff_mmap, threads, verbose)?;
    let t_compute = t3.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 4: Compute breadth: {:.2?}", t_compute);
    }

    // Step 5: write results
    let t4 = Instant::now();
    
    let out: Box<dyn Write> = match &args.output {
        Some(path) => {
            let file = File::create(path)?;
            Box::new(BufWriter::with_capacity(WRITE_BUF_SIZE, file))
        }
        None => {
            let stdout = std::io::stdout();
            let handle = stdout.lock();
            Box::new(BufWriter::with_capacity(WRITE_BUF_SIZE, handle))
        }
    };
    write_breadth_results(id_map, out, verbose)?;
    let t_write_out = t4.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 5: Write output: {:.2?}", t_write_out);
        let total = t0.elapsed();
        eprintln!("[TIMER] [run] Total time: {:.2?}", total);
    }

    Ok(())
}
