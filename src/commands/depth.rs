use anyhow::{Result, bail, Context};
use rayon::prelude::*;
use memmap2::Mmap;
use rust_htslib::bam::{self, Read};
use rust_htslib::bam::ext::BamRecordExtensions;
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

// Sentinel for missing entries
const MISSING: u64 = u64::MAX; 
// Number of IoSlices per batch writer
const WRITE_BUF_SIZE: usize = 32 * 1024 * 1024;
// BufWriter buffer size
const BATCH_SIZE: usize = 100_000;

/// Compute bin index for coordinate `x` given shift k.
/// Each bin has width 2^k bp. Smaller k = finer bins, larger k = coarser bins.
#[inline]
fn bin_of(x: u32, shift: u32) -> u32 {
    x >> shift
}

/// Arguments for `depth` command
#[derive(Parser, Debug)]
#[command(
    about = "Compute coverage depth across genomic features",
    long_about = "This tool computes sequencing depth (number of overlapping regions/reads per feature) \
                  from SAM/BAM/CRAM or BED input. It does not compute breadth/fraction coverage."
)]
pub struct DepthArgs {
    /// Input GFF file path
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Input source (BAM/SAM/CRAM or BED)
    #[arg(short = 's', long)]
    pub source: PathBuf,
    
    /// Output file (stdout if not provided)
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Bin width parameter (2^k bp) for spatial bucketing of features and queries.
    /// Choose k so that a typical read and feature span ~1–2 bins.
    ///
    /// Typical values:
    ///   Short reads (Illumina 100–150 bp): k=10–11 (1–2 kb bins)
    ///   PacBio HiFi (15–20 kb): k=13–14 (8–16 kb bins)
    ///   ONT long reads (30–60 kb): k=14–15 (16–32 kb bins)
    ///
    /// Adjust:
    ///   Longer features → increase k (larger bins, less fragmentation).
    ///   Denser features → decrease k (smaller bins, stronger filtering).
    #[arg(long = "bin-shift", default_value_t = 12)]
    pub bin_shift: u32,
    
    /// Number of threads for parallel processing
    #[arg(short = 't', long = "threads", default_value_t = 12)]
    pub threads: usize,

    /// Enable verbose output
    #[arg(short = 'v', long = "verbose", default_value_t = false)]
    pub verbose: bool,
}

/// Check half-open overlap: [a1, a2) vs [b1, b2)
#[inline(always)]
fn overlaps(a1: u32, a2: u32, b1: u32, b2: u32) -> bool {
    let left  = if a1 > b1 { a1 } else { b1 };
    let right = if a2 < b2 { a2 } else { b2 };
    left < right
}

/// Parse unsigned int quickly into u32
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

#[derive(Clone, Copy)]
struct RegionRef { start: u32, end: u32 }

#[derive(Clone, Copy)]
struct FeatureInst { start: u32, end: u32, id_idx: u32 }

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

/// Parse one GFF slice and count feature *depth* (how many regions overlap it; deduped per region)
fn compute_root_depth(
    gff_slice: &[u8],
    regions: &[RegionRef],
    bin_shift: u32,  
) -> FxHashMap<String, (String, u32, u32, usize)> {
    let mut id_to_idx: FxHashMap<String, u32> = FxHashMap::default();
    let mut id_strings: Vec<String> = Vec::new();
    let mut id_chrom: Vec<String> = Vec::new();
    let mut feats: Vec<FeatureInst> = Vec::new();

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
            let fstart0 = s1.saturating_sub(1);
            let fend0   = e1;

            if let Some(id) = fast_id(attrs) {
                let idx = *id_to_idx.entry(id.to_owned()).or_insert_with(|| {
                    let k = id_strings.len() as u32;
                    id_strings.push(id.to_owned());
                    id_chrom.push(seqid.to_owned());
                    k
                });
                feats.push(FeatureInst { start: fstart0, end: fend0, id_idx: idx });
            }
        }
    }

    if feats.is_empty() || regions.is_empty() {
        return FxHashMap::default();
    }

    let max_feat_bin = feats.iter().map(|f| bin_of(f.end.saturating_sub(1), bin_shift)).max().unwrap_or(0);
    let mut feat_bins: Vec<Vec<u32>> = vec![Vec::new(); (max_feat_bin as usize) + 1];
    for (i, f) in feats.iter().enumerate() {
        let b0 = bin_of(f.start, bin_shift);
        let b1 = bin_of(f.end.saturating_sub(1), bin_shift);
        for b in b0..=b1 {
            feat_bins[b as usize].push(i as u32);
        }
    }

    let n_ids = id_strings.len();
    let mut min_s: Vec<u32> = vec![u32::MAX; n_ids];
    let mut max_e: Vec<u32> = vec![0; n_ids];
    let mut depths: Vec<usize> = vec![0; n_ids];

    let mut cand: Vec<u32> = Vec::new();
    let mut hit_ids: Vec<u32> = Vec::new();

    for r in regions {
        cand.clear();
        let rb0 = bin_of(r.start, bin_shift);
        let rb1 = bin_of(r.end.saturating_sub(1), bin_shift);
        for b in rb0..=rb1 {
            if let Some(v) = feat_bins.get(b as usize) {
                cand.extend_from_slice(v);
            }
        }
        cand.sort_unstable();
        cand.dedup();

        hit_ids.clear();
        for &fi in &cand {
            let f = feats[fi as usize];
            if overlaps(f.start, f.end, r.start, r.end) {
                hit_ids.push(f.id_idx);
                if f.start < min_s[f.id_idx as usize] { min_s[f.id_idx as usize] = f.start; }
                if f.end   > max_e[f.id_idx as usize] { max_e[f.id_idx as usize] = f.end; }
            }
        }
        hit_ids.sort_unstable();
        hit_ids.dedup();
        for &ii in &hit_ids {
            depths[ii as usize] += 1;
        }
    }

    let mut out: FxHashMap<String, (String, u32, u32, usize)> = FxHashMap::default();
    for (i, &d) in depths.iter().enumerate() {
        if d > 0 {
            let s = if min_s[i] == u32::MAX { 0 } else { min_s[i] };
            out.insert(id_strings[i].clone(), (id_chrom[i].clone(), s, max_e[i], d));
        }
    }
    out
}

/// Batch API: for a batch of regions, return "feature ID -> (chrom, start, end, depth)".
///
/// - depth  = how many regions overlap with the feature (count of regions, deduped per region)
pub fn compute_hit_depth(
    index_data: &TreeIndexData,
    regions: &[(u32, u32, u32)],
    gof: &GofMap,
    gff_mmap: &Mmap,
    bin_shift: u32,
    threads: usize,
) -> Result<FxHashMap<String, (String, u32, u32, usize)>> {
    let mut by_root: FxHashMap<u32, Vec<RegionRef>> = FxHashMap::default();
    let idx = gof.index_cached();
    let gff_bytes: &[u8] = &gff_mmap[..];

    let mut hits: Vec<&Interval<u32>> = Vec::new();
    for &(chr, rstart, rend) in regions {
        if let Some(tree) = index_data.chr_entries.get(&chr) {
            hits.clear();
            tree.query_interval(rstart, rend, &mut hits);
            let mut seen_in_region: FxHashSet<u32> = FxHashSet::default();
            for h in &hits {
                if !seen_in_region.insert(h.root_fid) { continue; }
                if let Some(&(s_off, e_off)) = idx.get(&h.root_fid) {
                    if s_off == MISSING || e_off == MISSING || e_off <= s_off { continue; }
                    by_root.entry(h.root_fid).or_default().push(RegionRef { start: rstart, end: rend });
                }
            }
        }
    }

    let mut out: FxHashMap<String, (String, u32, u32, usize)> = FxHashMap::default();
    if by_root.is_empty() { return Ok(out); }

    let roots_iter = by_root.into_iter();
    if threads > 1 {
        // Parallel execution: each root slice is processed independently
        let partials: Vec<_> = roots_iter.par_bridge().map(|(root, regs)| {
            let (s_off, e_off) = *idx.get(&root).unwrap();
            let su = usize::try_from(s_off).unwrap();
            let eu = usize::try_from(e_off).unwrap();
            compute_root_depth(&gff_bytes[su..eu], &regs, bin_shift)
        }).collect();
        
        // Merge results from all roots
        for m in partials {
            for (id, (chrom, s, e, d)) in m {
                out.entry(id).and_modify(|(c0, s0, e0, depth)| {
                    if s < *s0 { *s0 = s; }
                    if e > *e0 { *e0 = e; }
                    *depth += d;
                    let _ = c0;
                }).or_insert((chrom, s, e, d));
            }
        }
    } else {
        // Serial execution
        for (root, regs) in roots_iter {
            let (s_off, e_off) = *idx.get(&root).unwrap();
            let su = usize::try_from(s_off).unwrap();
            let eu = usize::try_from(e_off).unwrap();
            let m = compute_root_depth(&gff_bytes[su..eu], &regs, bin_shift);
            for (id, (chrom, s, e, d)) in m {
                out.entry(id).and_modify(|(c0, s0, e0, depth)| {
                    if s < *s0 { *s0 = s; }
                    if e > *e0 { *e0 = e; }
                    *depth += d;
                    let _ = c0;
                }).or_insert((chrom, s, e, d));
            }
        }
    }

    Ok(out)
}

/// Process BAM input with mmap/htslib and batch queries.
/// Returns: "feature ID -> (chrom, start, end, depth)".
pub fn process_bam(
    bam_path: &Path,
    index_data: &TreeIndexData,
    gof: GofMap,
    gff_mmap: Mmap,
    bin_shift: u32,
    threads: usize,
    verbose: bool,
) -> Result<FxHashMap<String, (String, u32, u32, usize)>> {
    let mut global_id_counts: FxHashMap<String, (String, u32, u32, usize)> = FxHashMap::default();

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

    let mut batch: Vec<(u32, u32, u32)> = Vec::with_capacity(BATCH_SIZE);

    // Timers
    let mut t_parse = Duration::ZERO;
    let mut t_tidmap = Duration::ZERO;
    let mut t_filtermap = Duration::ZERO;
    let mut t_tree = Duration::ZERO;
    let mut t_depthmap = Duration::ZERO;

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
            if start_i64 < 0 {
                t_tidmap += t1.elapsed();
                continue;
            }
            let end_i64 = rec.reference_end();
            if end_i64 <= start_i64 {
                t_tidmap += t1.elapsed();
                continue;
            }
            let start = (start_i64 as i128).clamp(0, u32::MAX as i128) as u32;
            let end   = (end_i64 as i128).clamp(0, u32::MAX as i128) as u32;

            batch.push((chr_id, start, end));
        }
        t_tidmap += t1.elapsed();

        if batch.len() >= BATCH_SIZE {
            let t2 = Instant::now();
            let regions = std::mem::take(&mut batch);
            t_filtermap += t2.elapsed();

            let t3 = Instant::now();
            let id_counts = compute_hit_depth(index_data, &regions, &gof, &gff_mmap, bin_shift, threads)?;
            t_tree += t3.elapsed();

            let t4 = Instant::now();
            for (id, (chrom, s, e, d)) in id_counts {
                global_id_counts.entry(id).and_modify(|(c0, s0, e0, depth)| {
                    if s < *s0 { *s0 = s; }
                    if e > *e0 { *e0 = e; }
                    *depth += d;
                    let _ = c0;
                }).or_insert((chrom, s, e, d));
            }
            t_depthmap += t4.elapsed();
            batch.clear();
        }
    }

    if !batch.is_empty() {
        let t2 = Instant::now();
        let regions = std::mem::take(&mut batch);
        t_filtermap += t2.elapsed();

        let t3 = Instant::now();
        let id_counts = compute_hit_depth(index_data, &regions, &gof, &gff_mmap, bin_shift, threads)?;
        t_tree += t3.elapsed();

        let t4 = Instant::now();
        for (id, (chrom, s, e, d)) in id_counts {
            global_id_counts.entry(id).and_modify(|(c0, s0, e0, depth)| {
                if s < *s0 { *s0 = s; }
                if e > *e0 { *e0 = e; }
                *depth += d;
                let _ = c0;
            }).or_insert((chrom, s, e, d));
        }
        t_depthmap += t4.elapsed();
    }

    if verbose {
        eprintln!("[TIMER] (1) Opening BAM file:     {:.2?}", t_open_elapsed);
        eprintln!("[TIMER] (1b) Build tid2num map:  {:.2?}", t_map_build_elapsed);
        eprintln!("[TIMER] (2) Parsing records:     {:.2?}", t_parse);
        eprintln!("[TIMER] (3) Chrom ID mapping:    {:.2?}", t_tidmap);
        eprintln!("[TIMER] (4) Batch filter_map:    {:.2?}", t_filtermap);
        eprintln!("[TIMER] (5) Interval tree query: {:.2?}", t_tree);
        eprintln!("[TIMER] (6) DepthMap updates:    {:.2?}", t_depthmap);
    }

    Ok(global_id_counts)
}

/// Process BED input with mmap, parallel line parsing, and batch queries.
/// 
/// Returns: "feature ID -> (chrom, start, end, depth)".
/// - depth   = number of regions overlapping the feature
pub fn process_bed(
    bed_path: &Path,
    index_data: &TreeIndexData,
    gof: GofMap,
    gff_mmap: Mmap,
    bin_shift: u32,
    threads: usize,
    verbose: bool,
) -> Result<FxHashMap<String, (String, u32, u32, usize)>> {
    let mut global_id_counts: FxHashMap<String, (String, u32, u32, usize)> = FxHashMap::default();

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

    // process chunks in batches
    for chunk in line_offsets.windows(2).collect::<Vec<_>>().chunks(BATCH_SIZE) {
        // parse BED lines into (chr_id, start, end)
        let regions: Vec<(u32, u32, u32)> = chunk
            .par_iter()
            .filter_map(|w| {
                let start = w[0];
                let end = w[1];
                if start >= end {
                    return None;
                }
                let line = &data[start..end];
                if line.is_empty() || line[0] == b'#' {
                    return None;
                }

                let fields: Vec<&[u8]> = line
                    .split(|&b| b == b'\t' || b == b' ')
                    .filter(|f| !f.is_empty())
                    .collect();
                if fields.len() < 3 {
                    return None;
                }

                let chrom = std::str::from_utf8(fields[0]).ok()?;
                let s = std::str::from_utf8(fields[1]).ok()?.parse::<u32>().ok()?;
                let e = std::str::from_utf8(fields[2]).ok()?.trim_end().parse::<u32>().ok()?;
                if s >= e {
                    return None;
                }

                let &chr_num = index_data.seqid_to_num.get(chrom)?;
                Some((chr_num, s, e))
            })
            .collect();

        // compute depth only
        let id_counts = compute_hit_depth(index_data, &regions, &gof, &gff_mmap, bin_shift, threads)?;

        // merge into global results
        for (id, (chrom, s, e, d)) in id_counts {
            global_id_counts.entry(id).and_modify(|(c0, s0, e0, depth)| {
                if s < *s0 { *s0 = s; }
                if e > *e0 { *e0 = e; }
                *depth += d;
                let _ = c0; // chrom assumed consistent
            }).or_insert((chrom, s, e, d));
        }
    }

    Ok(global_id_counts)
}

/// Write "id\tchr\tstart\tend\tdepth" per line to output file.
pub fn write_depth_results<W: Write>(
    id_counts: FxHashMap<String, (String, u32, u32, usize)>,
    mut out: W,
    verbose: bool,
) -> Result<()> {
    use std::fmt::Write as FmtWrite;
    let mut buf = String::with_capacity(WRITE_BUF_SIZE);
    let mut written = 0usize;

    writeln!(buf, "id\tchr\tstart\tend\tdepth")?;
    
    for (id, (chr, start, end, depth)) in id_counts {
        writeln!(buf, "{id}\t{chr}\t{start}\t{end}\t{depth}")?;
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
        eprintln!("[INFO] Wrote {written} ID depth records");
    }
    Ok(())
}

/// Main entry for depth pipeline
pub fn run(args: &DepthArgs) -> Result<()> {
    let verbose = args.verbose;
    let bin_shift = args.bin_shift;
    let threads = if args.threads == 0 {
            std::thread::available_parallelism()
                .map(|n| n.get())
                .unwrap_or(1)
        } else {
            args.threads
        };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global();
    let gff_path = &args.input;
    
    // Step 1: load GFF index
    let t0 = Instant::now();
    let gof = load_gof(&gff_path)?;
    let file = File::open(gff_path).with_context(|| format!("Cannot open GFF file: {:?}", gff_path))?;
    let gff_mmap = unsafe { Mmap::map(&file) }.with_context(|| format!("GFF mmap failed for {:?}", gff_path))?;
    let t_load_fts = t0.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 1: Loading index took {:.2?}", t_load_fts);
    }

    // Step 2: build interval index
    let t1 = Instant::now();
    let index_data = TreeIndexData::load_tree_index(&gff_path)?;
    let t_build_index = t1.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 2: Building tree index took {:.2?}", t_build_index);
    }

    // Step 3: process input file
    let t2 = Instant::now();
    let source_path = &args.source;

    let ext = source_path
        .extension()
        .and_then(|s| s.to_str())
        .map(|s| s.to_lowercase());

    let id_counts = match ext.as_deref() {
        Some("bam") | Some("sam") | Some("cram") => {
            process_bam(source_path.as_path(), &index_data, gof, gff_mmap, bin_shift, threads, verbose)?
        }
        Some("bed") => {
            process_bed(source_path.as_path(), &index_data, gof, gff_mmap, bin_shift, threads, verbose)?
        }
        _ => {
            bail!(
                "Unsupported file type: {:?}. Expected .bam/.sam/.cram or .bed",
                source_path
            );
        }
    };
    let t_process_input = t2.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 3: Processing input took {:.2?}", t_process_input);
    }

    // Step 4: write results
    let t3 = Instant::now();
    
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
    
    write_depth_results(id_counts, out, verbose)?;
    
    let t_write_out = t3.elapsed();
    if verbose {
        eprintln!("[TIMER] [run] Step 4: Writing results took {:.2?}", t_write_out);
    }

    if verbose {
        let total = t0.elapsed();
        eprintln!("[TIMER] [run] Total pipeline time: {:.2?}", total);
    }

    Ok(())
}