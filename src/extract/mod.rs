use crate::{
    CommonArgs, extract_root_matches, load_fts, load_gof, load_prt, roots_to_offsets,
    write_gff_output,
};
use anyhow::{Context, Result, bail};
use clap::Parser;
use memchr::{memchr, memmem};
use memmap2::Mmap;
use rustc_hash::{FxHashMap, FxHashSet};
use std::{
    fs::File,
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
    //   sync::{Arc, Mutex},
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

/// When `--match-only` is ON: within each root block, emit only lines whose `ID` exactly matches
/// the user-specified features under that root. Optional `types_filter` is applied to column 3.
/// - `fts` is used to map numeric feature IDs back to string IDs for matching.
fn write_gff_output_filtered(
    gff_path: &PathBuf,
    blocks: &[(u32, u64, u64)],
    per_root_matches: &FxHashMap<u32, FxHashSet<u32>>,
    fts: &[String],
    output_path: &Option<PathBuf>,
    types_filter: Option<&str>,
    verbose: bool,
) -> Result<()> {
    let file =
        File::open(gff_path).with_context(|| format!("Cannot open GFF file: {:?}", gff_path))?;
    let mmap =
        unsafe { Mmap::map(&file) }.with_context(|| format!("mmap failed for {:?}", gff_path))?;
    let file_len = mmap.len();

    let type_allow: Option<FxHashSet<String>> = types_filter.map(|s| {
        s.split(',')
            .map(|t| t.trim().to_string())
            .filter(|t| !t.is_empty())
            .collect()
    });

    use rayon::prelude::*;
    let mut parts: Vec<(u64, Vec<u8>)> = blocks
        .par_iter()
        .filter_map(|&(root, start, end)| {
            let matched_numeric = per_root_matches.get(&root)?;
            if matched_numeric.is_empty() {
                return None;
            }
            let keep: FxHashSet<&str> = matched_numeric
                .iter()
                .filter_map(|&n| fts.get(n as usize))
                .map(|s| s.as_str())
                .collect();
            if keep.is_empty() {
                return None;
            }

            let s = start as usize;
            let e = end.min(file_len as u64) as usize;
            if s >= e || e > file_len {
                return None;
            }
            let window = &mmap[s..e];

            let mut out = Vec::<u8>::with_capacity(1024);
            let mut pos = 0usize;

            let next_line = |from: usize| -> Option<(usize, usize, usize)> {
                if from >= window.len() {
                    return None;
                }
                let rel = memchr(b'\n', &window[from..])
                    .map(|i| from + i + 1)
                    .unwrap_or(window.len());
                let mut end_no_nl = rel;
                if end_no_nl > from && window[end_no_nl - 1] == b'\n' {
                    end_no_nl -= 1;
                }
                if end_no_nl > from && window[end_no_nl - 1] == b'\r' {
                    end_no_nl -= 1;
                }
                Some((from, rel, end_no_nl))
            };

            let type_ok = |line: &[u8]| -> bool {
                if let Some(allow) = &type_allow {
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
                        return allow.contains(ty_str);
                    } else {
                        return false;
                    }
                }
                true
            };

            let id_hits_keep = |line_no_crlf: &[u8]| -> bool {
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
                if let Some(p) = memmem::find(attr, b"ID=") {
                    let vstart = p + 3;
                    // ID 值以 ';' 或行尾结束
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

            while let Some((ls, le, ln_end)) = next_line(pos) {
                pos = le;
                let line = &window[ls..le];
                if !line.is_empty() && line[0] == b'#' {
                    continue;
                }
                let line_no_crlf = &window[ls..ln_end];

                if !type_ok(line_no_crlf) {
                    continue;
                }
                // ID 命中
                if id_hits_keep(line_no_crlf) {
                    out.extend_from_slice(line);
                }
            }

            if verbose {
                eprintln!(
                    "[filter] root={} block=[{}..{}] keep_ids={} matched_lines={}",
                    root,
                    start,
                    end,
                    keep.len(),
                    out.iter().filter(|&&b| b == b'\n').count()
                );
            }

            if out.is_empty() {
                None
            } else {
                Some((start, out))
            }
        })
        .collect();

    parts.sort_unstable_by_key(|(s, _)| *s);

    let raw: Box<dyn Write> = match output_path {
        Some(p) => {
            Box::new(File::create(p).with_context(|| format!("Cannot create output: {:?}", p))?)
        }
        None => Box::new(std::io::stdout()),
    };
    let mut writer = BufWriter::new(raw);
    for (_, buf) in parts {
        writer.write_all(&buf)?;
    }
    writer.flush()?;

    Ok(())
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

    // Build lookup maps
    let fts_map: FxHashMap<&str, u32> = fts
        .iter()
        .enumerate()
        .map(|(i, s)| (s.as_str(), i as u32))
        .collect();
    let prt_map: FxHashMap<u32, u32> = prt.iter().map(|e| (e.child, e.parent)).collect();
    let gof_map: FxHashMap<u32, (u64, u64)> = gof
        .iter()
        .map(|e| (e.feature_id, (e.start_offset, e.end_offset)))
        .collect();

    // Read feature IDs
    let feature_ids: FxHashSet<String> = if let Some(ref file_path) = args.feature_file {
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
    let root_matches = extract_root_matches(
        &feature_ids,
        &fts_map,
        &prt_map,
        args.common.effective_threads(),
    );

    // Extract unique roots, then Phase B: map to offsets
    let roots: Vec<u32> = {
        let mut v: Vec<u32> = root_matches.iter().map(|rm| rm.root).collect();
        v.sort_unstable();
        v.dedup();
        v
    };
    let root_blocks: Vec<(u32, u64, u64)> = roots_to_offsets(&roots, &gof_map);

    if !args.common.full_model {
        // Build per-root matched numeric ID sets
        let per_root_matches: FxHashMap<u32, FxHashSet<u32>> = {
            let mut m: FxHashMap<u32, FxHashSet<u32>> = FxHashMap::default();
            m.reserve(root_matches.len());
            for rm in &root_matches {
                m.entry(rm.root)
                    .or_default()
                    .extend(rm.matched.iter().copied());
            }
            m
        };

        // NEW: only emit exactly matched lines within blocks
        write_gff_output_filtered(
            gff_path,
            &root_blocks,
            &per_root_matches,
            &fts,
            &args.common.output,
            args.common.types.as_deref(),
            verbose,
        )?;
    } else {
        // Back-compat path: emit full blocks using the existing function
        let block_offsets: Vec<(u64, u64)> = root_blocks.iter().map(|&(_, s, e)| (s, e)).collect();

        write_gff_output(
            gff_path,
            &block_offsets,
            &args.common.output,
            args.common.types.as_deref(),
            verbose,
        )?;
    }

    if verbose {
        eprintln!("[timing] Total elapsed: {:?}", overall_start.elapsed());
    }

    Ok(())
}
