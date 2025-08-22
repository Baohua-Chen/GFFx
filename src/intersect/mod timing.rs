use anyhow::{Context, Result};
use clap::{ArgGroup, Parser};
use memchr::memchr;
use memmap2::Mmap;
use rayon::prelude::*;
use std::{
    fs::File,
    io::{self, BufReader, BufWriter, IoSlice, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
    time::Instant,
};
//use byteorder::{LittleEndian, ReadBytesExt};
use bincode2::deserialize_from;
use lexical_core::parse;
use rustc_hash::{FxHashMap, FxHashSet};
//use serde::{Serialize, de::DeserializeOwned};

use crate::{
    CommonArgs, Interval, IntervalTree, Timing, append_suffix, load_gof, load_sqs,
    roots_to_offsets, write_gff_output,
};

/// 批量 writev 的 IoSlice 数量（可按机器微调做基准）
const IOV_BATCH: usize = 256;
/// BufWriter 缓冲大小（可按机器微调做基准）
const WRITE_BUF_SIZE: usize = 32 * 1024 * 1024;

/// 补齐：命中的 root 及其下属 feature（matched 先占位）
#[derive(Debug, Clone)]
pub struct RootMatched {
    pub root: u32,
    pub matched: Vec<u32>,
}

/// Arguments for region intersection operations
#[derive(Parser, Debug)]
#[command(
    about = "Extract models by a region or regions from a BED file",
    long_about = "This tool extracts features and their parent models that intersect with specified regions"
)]
#[clap(group(
    ArgGroup::new("regions").required(true).args(&["region", "bed"])
))]
#[clap(group(
    ArgGroup::new("mode").args(&["contained", "contains_region", "overlap"])
))]
pub struct IntersectArgs {
    #[clap(flatten)]
    pub common: CommonArgs,

    /// Single region in format "chr:start-end"
    #[arg(short = 'r', long, group = "regions")]
    pub region: Option<String>,

    /// BED file containing regions
    #[arg(short = 'b', long, group = "regions")]
    pub bed: Option<PathBuf>,

    /// Only return features fully contained within regions
    #[arg(short = 'c', long, group = "mode")]
    pub contained: bool,

    /// Only return features that fully contain the regions
    #[arg(short = 'C', long, group = "mode")]
    pub contains_region: bool,

    /// Return any overlapping features (default)
    #[arg(short = 'O', long, group = "mode")]
    pub overlap: bool,

    /// Invert the selection (exclude matching features)
    #[arg(short = 'I', long, default_value_t = false)]
    pub invert: bool,

    #[arg(short = 'M', long = "match-only", default_value_t = false)]
    pub match_only: bool,
}

/// Index data structure using interval trees
#[derive(Debug)]
pub struct IndexData {
    pub chr_entries: FxHashMap<u32, IntervalTree<u32>>,
    pub seqid_to_num: FxHashMap<String, u32>,
    pub common: CommonArgs,
}

impl IndexData {
    pub fn load(index_prefix: &Path, common: &CommonArgs) -> Result<Self> {
        let (seqids, _) = load_sqs(index_prefix)?;
        let seqid_to_num = seqids
            .into_iter()
            .enumerate()
            .map(|(i, sq)| (sq, i as u32))
            .collect();

        let rit_path = append_suffix(index_prefix, ".rit");
        let rix_path = append_suffix(index_prefix, ".rix");
        let chr_entries = Self::load_region_index(&rit_path, &rix_path)?;

        Ok(Self {
            chr_entries,
            seqid_to_num,
            common: common.clone(),
        })
    }

    fn load_region_index(
        rit_path: &Path,
        rix_path: &Path,
    ) -> Result<FxHashMap<u32, IntervalTree<u32>>> {
        let mut reader = BufReader::new(File::open(rit_path)?);
        let offsets: Vec<u64> = {
            let f = File::open(rix_path)?;
            serde_json::from_reader(f)?
        };
        let mut map = FxHashMap::with_capacity_and_hasher(offsets.len(), Default::default());
        for (i, &offset) in offsets.iter().enumerate() {
            reader.seek(SeekFrom::Start(offset))?;
            let tree: IntervalTree<u32> = deserialize_from(&mut reader)?;
            map.insert(i as u32, tree);
        }
        Ok(map)
    }
}

/// Overlap detection modes
#[derive(Debug, Clone, Copy)]
enum OverlapMode {
    Contained,
    ContainsRegion,
    Overlap,
}

fn gff_type_allowed(line: &[u8], allow: &FxHashSet<String>) -> bool {
    // Fast parse the 3rd field (type) without allocations
    let mut off = 0usize;
    let mut tabs = 0u8;
    while tabs < 2 {
        match memchr(b'\t', &line[off..]) {
            Some(i) => {
                off += i + 1;
                tabs += 1;
            }
            None => return false,
        }
    }
    let i2 = match memchr(b'\t', &line[off..]) {
        Some(i) => off + i,
        None => return false,
    };
    let ty = &line[off..i2];
    match std::str::from_utf8(ty) {
        Ok(s) => allow.contains(s),
        Err(_) => false,
    }
}

/// Core feature query logic using interval trees
fn query_features(
    index_data: &IndexData,
    regions: Vec<(u32, u32, u32)>,
    contained: bool,
    contains_region: bool,
    invert: bool,
) -> Result<Vec<(u32, u32, u32)>> {
    let verbose = index_data.common.verbose;
    let t = Timing::new(verbose);

    if verbose {
        eprintln!(
            "[DEBUG] Starting query_features with {} regions",
            regions.len()
        );
        eprintln!(
            "[DEBUG] Mode: contained={}, contains_region={}, invert={}",
            contained, contains_region, invert
        );
    }

    let buckets: Vec<Vec<(u32, u32, u32)>> = {
        let _g = t.scoped("bucketize_regions");
        let mut b = vec![Vec::new(); index_data.seqid_to_num.len()];
        for (chr, start, end) in regions {
            b[chr as usize].push((chr, start, end));
        }
        b
    };

    let _mode = if contained {
        OverlapMode::Contained
    } else if contains_region {
        OverlapMode::ContainsRegion
    } else {
        OverlapMode::Overlap
    };

    let mut results = Vec::new();
    let mut total_hits: u128 = 0;

    {
        let _g = t.scoped("scan_interval_trees");
        for (&seq_num, tree) in &index_data.chr_entries {
            let chr_regs = &buckets[seq_num as usize];
            if chr_regs.is_empty() {
                continue;
            }
            if verbose {
                eprintln!(
                    "[DEBUG] Querying chromosome {} with {} regions",
                    seq_num,
                    chr_regs.len()
                );
            }
            for &(_, rstart, rend) in chr_regs {
                let hits: Vec<&Interval<u32>> = tree.query_interval(rstart, rend);
                total_hits += hits.len() as u128;
                results.extend(hits.into_iter().map(|iv| (iv.root_fid, iv.start, iv.end)));
            }
        }
    }

    t.set("total_hits", total_hits);
    t.set("chromosomes", index_data.chr_entries.len() as u128);
    t.finish("query_features");
    Ok(results)
}

/// Parse a single genomic region string (chr:start-end)
fn parse_region(
    region: &str,
    seqid_map: &FxHashMap<String, u32>,
    common: &CommonArgs,
) -> Result<(u32, u32, u32)> {
    let (seq, range) = region
        .split_once(':')
        .context("Invalid region format, expected 'chr:start-end'")?;
    let (s, e) = range
        .split_once('-')
        .context("Invalid range format, expected 'start-end'")?;
    let start = s.parse::<u32>()?;
    let end = e.parse::<u32>()?;
    let chr = seqid_map
        .get(seq)
        .with_context(|| format!("Sequence ID not found: {}", seq))?;
    if start >= end {
        anyhow::bail!("Region start must be less than end ({} >= {})", start, end);
    }
    if common.verbose {
        eprintln!(
            "[DEBUG] Parsed region: chr={}, start={}, end={}",
            chr, start, end
        );
    }
    Ok((*chr, start, end))
}

/// Parse BED file using mmap zero-copy field splitting
fn parse_bed_file(
    bed_path: &Path,
    seqid_map: &FxHashMap<String, u32>,
    common: &CommonArgs,
) -> Result<Vec<(u32, u32, u32)>> {
    let t = Timing::new(common.verbose);
    let mmap = {
        let _g = t.scoped("mmap_bed");
        let file = File::open(bed_path)?;
        unsafe { Mmap::map(&file)? }
    };
    let regions = {
        let _g = t.scoped("parse_lines");
        let mut regions = Vec::new();
        let mut lines: u128 = 0;
        for line in mmap.split(|&b| b == b'\n') {
            if line.is_empty() || line[0] == b'#' {
                continue;
            }
            lines += 1;
            let line_str = std::str::from_utf8(line)?;
            let mut parts = line_str.split_ascii_whitespace();
            let (Some(seq), Some(s), Some(e)) = (parts.next(), parts.next(), parts.next()) else {
                continue;
            };
            let Some(&chr) = seqid_map.get(seq) else {
                continue;
            };
            let start = parse::<u32>(s.as_bytes())?;
            let end = parse::<u32>(e.as_bytes())?;
            regions.push((chr, start, end));
        }
        t.set("bed_lines", lines);
        regions
    };
    if common.verbose {
        t.set("regions_parsed", regions.len() as u128);
        t.finish("parse_bed_file");
    }
    Ok(regions)
}

/// 仅按坐标过滤并输出（不依赖 ID），第三参是“查询区间列表映射”
///
/// 改造要点：
/// 1) Per-block parallel scan 收集 **(line_start, line_end)** 偏移（不复制字节）
/// 2) 对块按 `start` 排序；顺序遍历块，将行切片打包成 `IoSlice` 批次
/// 3) 用 `write_vectored` 批量写出（含部分写处理）
/// 4) 使用大号 `BufWriter`
fn write_gff_match_only_by_coords(
    gff_path: &Path,
    blocks: &[(u32, u64, u64)], // 每个 seqid -> 若干 [start,end] 闭区间
    query_ivmap: &FxHashMap<String, Vec<(u32, u32)>>,
    types_filter: Option<&str>,
    output_path: &Option<PathBuf>,
    verbose: bool,
) -> Result<()> {
    let t = Timing::new(verbose);

    // mmap the whole GFF once
    let (mmap, file_len) = {
        let _g = t.scoped("mmap_gff");
        let file = std::fs::File::open(gff_path)
            .with_context(|| format!("Cannot open GFF: {:?}", gff_path))?;
        let mmap = unsafe { Mmap::map(&file) }
            .with_context(|| format!("mmap failed for {:?}", gff_path))?;
        let len = mmap.len();
        (mmap, len)
    };

    // parse type filters to a set once
    let type_allow: Option<FxHashSet<String>> = {
        let _g = t.scoped("parse_types_filter");
        types_filter.map(|s| {
            s.split(',')
                .map(|t| t.trim().to_string())
                .filter(|t| !t.is_empty())
                .collect()
        })
    };

    // Parallel scan blocks: produce (block_start, Vec<(line_start,line_end)>)
    // Note: we never copy line bytes, only collect offsets.
    let mut parts: Vec<(u64, Vec<(u64, u64)>)> = {
        let _g = t.scoped("scan_blocks_parallel");
        let bytes_out = std::sync::atomic::AtomicU64::new(0);

        let parts: Vec<(u64, Vec<(u64, u64)>)> = blocks
            .par_iter()
            .filter_map(|&(_root, start, end)| {
                let s = start as usize;
                let e = (end as usize).min(file_len);
                if s >= e || e > file_len {
                    return None;
                }
                let src = &mmap[s..e];

                // Collect matched line ranges as global file offsets.
                let mut matched_offsets: Vec<(u64, u64)> = Vec::with_capacity(256);
                let mut pos = 0usize;

                while pos < src.len() {
                    // Find next newline boundary
                    let nl = match memchr(b'\n', &src[pos..]) {
                        Some(i) => pos + i + 1, // include '\n'
                        None => src.len(),
                    };
                    let line = &src[pos..nl];

                    // Trim trailing '\n' for parsing
                    let line_nocr = if line.ends_with(b"\n") {
                        &line[..line.len() - 1]
                    } else {
                        line
                    };

                    if !line_nocr.is_empty() && line_nocr[0] != b'#' {
                        // Optional: type filter first to early discard
                        let mut pass = true;
                        if let Some(allow) = &type_allow
                            && !gff_type_allowed(line_nocr, allow)
                        {
                            pass = false;
                        }
                        if pass && gff_line_overlaps_queries(line_nocr, query_ivmap) {
                            // Record absolute offsets in the file (including '\n')
                            let abs_start = start + pos as u64;
                            let abs_end = start + nl as u64;
                            // Safety: bounds already clamped by file_len
                            matched_offsets.push((abs_start, abs_end));
                            bytes_out.fetch_add(
                                (abs_end - abs_start) as u64,
                                std::sync::atomic::Ordering::Relaxed,
                            );
                        }
                    }

                    pos = nl;
                }

                if matched_offsets.is_empty() {
                    None
                } else {
                    Some((start, matched_offsets))
                }
            })
            .collect();

        t.set(
            "bytes_selected",
            bytes_out.load(std::sync::atomic::Ordering::Relaxed) as u128,
        );
        t.set("blocks_in", blocks.len() as u128);
        parts
    };

    // Keep global order stable by block start (we don't merge offsets as per user's requirement)
    {
        let _g = t.scoped("sort_parts");
        parts.sort_unstable_by_key(|(s, _)| *s);
    }

    // Helper: write all slices using write_vectored with partial-write handling.
    // We construct a temporary Vec<IoSlice> per batch; batch size is small (<= IOV_BATCH).
    fn write_all_vectored<W: Write>(w: &mut W, mut slices: Vec<&[u8]>) -> io::Result<()> {
        // Fast path: nothing to write
        if slices.is_empty() {
            return Ok(());
        }

        // Keep writing until all slices are fully consumed
        while !slices.is_empty() {
            // Rebuild IoSlice views for current remainder
            let iov: Vec<IoSlice<'_>> = slices.iter().map(|s| IoSlice::new(s)).collect();

            let wrote = w.write_vectored(&iov)?;
            if wrote == 0 {
                return Err(io::Error::new(
                    io::ErrorKind::WriteZero,
                    "write_vectored returned 0",
                ));
            }

            // Consume 'wrote' bytes from the front of `slices`
            let mut remaining = wrote;
            let mut drop_count = 0;

            for s in &mut slices {
                if remaining == 0 {
                    break;
                }
                if remaining >= s.len() {
                    remaining -= s.len();
                    drop_count += 1;
                } else {
                    // Advance within the first partially-written slice
                    *s = &s[remaining..];
                    remaining = 0;
                }
            }

            if drop_count > 0 {
                slices.drain(0..drop_count);
            }
        }

        Ok(())
    }

    // Write out: use large BufWriter and batch IoSlice slices across consecutive parts.
    {
        let _g = t.scoped("write_output");

        // Assemble and write batches, reusing a small Vec<&[u8]> to avoid reallocs
        let mut batch: Vec<&[u8]> = Vec::with_capacity(IOV_BATCH);

        if let Some(p) = output_path {
            // File output path: create file and large BufWriter
            let file = std::fs::File::create(p)?;
            let mut writer = BufWriter::with_capacity(WRITE_BUF_SIZE, file);

            for (_, ranges) in parts.iter() {
                for &(ls, le) in ranges {
                    // Safety: ls/le were validated against file_len earlier
                    let slice = &mmap[ls as usize..le as usize];
                    batch.push(slice);
                    if batch.len() >= IOV_BATCH {
                        write_all_vectored(&mut writer, std::mem::take(&mut batch))?;
                    }
                }
            }
            if !batch.is_empty() {
                write_all_vectored(&mut writer, std::mem::take(&mut batch))?;
            }
            writer.flush()?;
        } else {
            // Stdout path: lock stdout and use large BufWriter
            let stdout = std::io::stdout();
            let handle = stdout.lock();
            let mut writer = BufWriter::with_capacity(WRITE_BUF_SIZE, handle);

            for (_, ranges) in parts.iter() {
                for &(ls, le) in ranges {
                    let slice = &mmap[ls as usize..le as usize];
                    batch.push(slice);
                    if batch.len() >= IOV_BATCH {
                        write_all_vectored(&mut writer, std::mem::take(&mut batch))?;
                    }
                }
            }
            if !batch.is_empty() {
                write_all_vectored(&mut writer, std::mem::take(&mut batch))?;
            }
            writer.flush()?;
        }
    }

    if verbose {
        t.finish("");
    }
    if verbose {
        eprintln!(
            "[INFO] match-only by coords 完成；输入分块 {}",
            blocks.len()
        );
    }
    Ok(())
}

/// 解析 GFF 的第 3 列 type，并做集合判断。
fn gff_line_overlaps_queries(line: &[u8], ivmap: &FxHashMap<String, Vec<(u32, u32)>>) -> bool {
    // Parse first 5 columns quickly: seq, source, type, start, end
    let mut off = 0usize;

    let i1 = match memchr(b'\t', &line[off..]) {
        Some(i) => off + i,
        None => return false,
    };
    let seq = &line[off..i1];
    off = i1 + 1;

    let i2 = match memchr(b'\t', &line[off..]) {
        Some(i) => off + i,
        None => return false,
    };
    off = i2 + 1;

    let i3 = match memchr(b'\t', &line[off..]) {
        Some(i) => off + i,
        None => return false,
    };
    off = i3 + 1;

    let i4 = match memchr(b'\t', &line[off..]) {
        Some(i) => off + i,
        None => return false,
    };
    let start = match parse_u32_ascii(&line[off..i4]) {
        Some(v) => v,
        None => return false,
    };
    off = i4 + 1;

    let i5 = match memchr(b'\t', &line[off..]) {
        Some(i) => off + i,
        None => return false,
    };
    let end = match parse_u32_ascii(&line[off..i5]) {
        Some(v) => v,
        None => return false,
    };

    let seq_str = match std::str::from_utf8(seq) {
        Ok(s) => s,
        Err(_) => return false,
    };
    let ivs = match ivmap.get(seq_str) {
        Some(v) => v,
        None => return false,
    };

    for &(qs, qe) in ivs {
        // Use simple integer comparisons to check overlap
        if (qs <= start && start <= qe)
            || (qs <= end && end <= qe)
            || (start <= qs && qs <= end)
            || (start <= qe && qe <= end)
        {
            return true;
        }
    }
    false
}

/// 零分配地从 ASCII 数字解析 u32
#[inline]
fn parse_u32_ascii(s: &[u8]) -> Option<u32> {
    let mut v: u32 = 0;
    if s.is_empty() {
        return None;
    }
    for &c in s {
        if !c.is_ascii_digit() {
            return None;
        }
        v = v.checked_mul(10)?.checked_add((c - b'0') as u32)?;
    }
    Some(v)
}

/// Main execution function
/// Main execution function
pub fn run(args: &IntersectArgs) -> Result<()> {
    let total = Instant::now();
    let verbose = args.common.verbose;
    let t = Timing::new(verbose);

    if verbose {
        eprintln!(
            "[INFO] Loading regions from: {}",
            args.common.input.display()
        );
    }

    let (_, seqid_map) = {
        let _g = t.scoped("load_sqs_prefix");
        load_sqs(&args.common.input)?
    };

    let regions = {
        let _g = t.scoped("parse_regions");
        if let Some(bed) = &args.bed {
            parse_bed_file(bed, &seqid_map, &args.common)?
        } else if let Some(r) = &args.region {
            vec![parse_region(r, &seqid_map, &args.common)?]
        } else {
            anyhow::bail!("No region specified");
        }
    };
    t.set("regions", regions.len() as u128);

    let index_data = {
        let _g = t.scoped("load_index_data");
        IndexData::load(&args.common.input, &args.common)?
    };

    let feats = {
        let _g = t.scoped("query_features");
        query_features(
            &index_data,
            regions.clone(),
            args.contained,
            args.contains_region,
            args.invert,
        )?
    };
    t.set("feature_hits", feats.len() as u128);

    // 收集 feature 的 ID（root_fid）
    let ids: Vec<u32> = {
        let _g = t.scoped("collect_dedup_ids");
        let mut ids: Vec<u32> = feats.iter().map(|&(id, _, _)| id).collect();
        ids.sort_unstable();
        ids.dedup();
        t.set("ids_unique", ids.len() as u128);
        ids
    };

    let gof_map: FxHashMap<_, _> = {
        let _g = t.scoped("load_gof_build_map");
        let gof = load_gof(&args.common.input)?;
        gof.into_iter()
            .map(|e| (e.feature_id, (e.start_offset, e.end_offset)))
            .collect()
    };

    if args.match_only {
        let blocks_and_write = Instant::now();

        let root_matches: Vec<RootMatched> = {
            let _g = t.scoped("group_roots_match_only");
            let mut grouped: FxHashMap<u32, Vec<u32>> = FxHashMap::default();
            for &(root, _s, _e) in &feats {
                grouped.entry(root).or_default().push(root);
            }
            grouped
                .into_iter()
                .map(|(root, matched)| RootMatched { root, matched })
                .collect()
        };
        t.set("root_groups", root_matches.len() as u128);

        let roots: Vec<u32> = {
            let _g = t.scoped("collect_roots");
            let mut s: FxHashSet<u32> = FxHashSet::default();
            for rm in &root_matches {
                s.insert(rm.root);
            }
            s.into_iter().collect()
        };

        let blocks: Vec<(u32, u64, u64)> = {
            let _g = t.scoped("roots_to_offsets");
            roots_to_offsets(&roots, &gof_map)
        };
        t.set("blocks", blocks.len() as u128);

        // Build query interval map by seq name
        let query_ivmap: FxHashMap<String, Vec<(u32, u32)>> = {
            let _g = t.scoped("build_query_ivmap");
            let mut num_to_seq: FxHashMap<u32, String> = FxHashMap::default();
            for (name, &num) in index_data.seqid_to_num.iter() {
                num_to_seq.insert(num, name.clone());
            }
            let mut m: FxHashMap<String, Vec<(u32, u32)>> = FxHashMap::default();
            for &(chr_num, s, e) in &regions {
                if let Some(seq_name) = num_to_seq.get(&chr_num) {
                    m.entry(seq_name.clone()).or_default().push((s, e));
                }
            }
            m
        };

        {
            let _g = t.scoped("write_match_only_by_coords");
            write_gff_match_only_by_coords(
                args.common.input.as_path(),
                &blocks,
                &query_ivmap,
                args.common.types.as_deref(),
                &args.common.output,
                args.common.verbose,
            )?;
        }

        if verbose {
            eprintln!(
                "[TIME] match_only pipeline: {:.3} ms",
                blocks_and_write.elapsed().as_secs_f64() * 1e3
            );
        }
    } else {
        // 非 -M 模式：沿用原有 write_gff_output
        let blocks: Vec<(u64, u64)> = {
            let _g = t.scoped("ids_to_blocks");
            let mut blocks: Vec<(u64, u64)> = Vec::with_capacity(ids.len());
            for &id in &ids {
                if let Some(&(s, e)) = gof_map.get(&id) {
                    if s < e {
                        blocks.push((s, e));
                    } else if verbose {
                        eprintln!("[WARNING] Invalid offsets for feature {}: {}..{}", id, s, e);
                    }
                } else if verbose {
                    eprintln!("[WARNING] Feature ID {} not found in GOF map", id);
                }
            }
            t.set("blocks_raw", blocks.len() as u128);
            blocks
        };

        {
            let _g = t.scoped("write_gff_output");
            write_gff_output(
                args.common.input.as_path(),
                &blocks,
                &args.common.output,
                args.common.types.as_deref(),
                args.common.verbose,
            )?;
        }
    }

    if verbose {
        t.finish("run_intersect");
        eprintln!("[TIME] Total: {:.2?}", total.elapsed());
    }
    Ok(())
}
