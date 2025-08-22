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

pub fn append_suffix(path: &Path, suffix: &str) -> PathBuf {
    let parent = path.parent().unwrap_or_else(|| Path::new(""));
    let filename = path.file_name().unwrap_or_default().to_string_lossy();
    parent.join(format!("{filename}{suffix}"))
}

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
        default_value_t = 8,
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
pub fn write_gff_output_old(
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

pub fn write_gff_output(
    gff_path: &Path,
    blocks: &[(u64, u64)],
    output_path: &Option<std::path::PathBuf>,
    _allowed_types: Option<&str>,
    verbose: bool,
) -> Result<()> {
    // 1) mmap
    let file = File::open(gff_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let file_len = mmap.len();

    // 2) sort + merge
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

    // 3) 组装切片（同时做边界校验）
    let mut slices: Vec<IoSlice<'_>> = Vec::with_capacity(merged.len());
    for &(so, eo) in &merged {
        // 边界检查
        if so >= eo {
            continue;
        }
        let (start, end) = (so as usize, eo as usize);
        if end > file_len {
            continue;
        }
        slices.push(IoSlice::new(&mmap[start..end]));
    }

    // 4) 打开输出
    let mut writer: Box<dyn Write> = match output_path {
        Some(p) => Box::new(BufWriter::new(File::create(p)?)),
        None => Box::new(BufWriter::new(stdout())),
    };

    // 5) 分批 + 处理部分写
    //    为了正确处理部分写，我们维护 (batch_index, slice_index, intra_offset)
    const MAX_IOV: usize = 1024; // 保守上限
    let mut base = 0;
    while base < slices.len() {
        let end = (base + MAX_IOV).min(slices.len());
        // 我们需要一个可变的“当前批次”指针序列，但 IoSlice 是不可变引用；
        // partial write 时对第一个未完全写完的切片，改用标量 write 写剩余部分，然后再进入下一切片。
        let batch = &slices[base..end];

        // 先尝试一次 vectored 写
        let nw = writer.write_vectored(batch)?;
        let mut remaining = nw;
        let mut i = 0;

        // 计算完整消费了多少个切片
        while i < batch.len() && remaining >= batch[i].len() {
            remaining -= batch[i].len();
            i += 1;
        }

        // 对第 i 个切片，如果有部分写，写完它的剩余部分
        if i < batch.len() && remaining > 0 {
            let cur = &batch[i];
            // 已写 remaining，写剩余 cur.len()-remaining
            writer.write_all(&cur[remaining..])?;
            i += 1;
        }

        // 继续把本批余下的切片逐个写完（不再用 vectored，避免再次处理复杂 partial write）
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

pub struct RootMatched {
    pub root: u32,
    pub matched: Vec<u32>,
}

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
