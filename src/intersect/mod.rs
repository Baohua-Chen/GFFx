use anyhow::{Result, Context};
use clap::{Parser, ArgGroup};
use memmap2::Mmap;
use std::{
    fs::File,
    io::{Write, Seek, SeekFrom, BufReader},
    path::{Path, PathBuf},
};
//use byteorder::{LittleEndian, ReadBytesExt};
use rustc_hash::FxHashMap;
use lexical_core::parse;
use bincode2::deserialize_from;
//use serde::{Serialize, de::DeserializeOwned};
use crate::{CommonArgs, load_sqs, append_suffix, load_gof, IntervalTree, Interval};

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
}

/// Overlap detection modes
#[derive(Debug, Clone, Copy)]
enum OverlapMode {
    Contained,
    ContainsRegion,
    Overlap,
}

/// Core feature query logic using interval trees
fn query_features(
    index_data: &IndexData,
    regions: Vec<(u32, u32, u32)>,
    contained: bool,
    contains_region: bool,
    invert: bool,
) -> Result<Vec<(u32, u32, u32)>> {
    if index_data.common.verbose {
        eprintln!("[DEBUG] Starting query_features with {} regions", regions.len());
        eprintln!("[DEBUG] Mode: contained={}, contains_region={}, invert={}",
                 contained, contains_region, invert);
    }

    let mut buckets: Vec<Vec<(u32, u32, u32)>> = vec![Vec::new(); index_data.seqid_to_num.len()];
    for (chr, start, end) in regions {
        buckets[chr as usize].push((chr, start, end));
    }

    let _mode = if contained {
        OverlapMode::Contained
    } else if contains_region {
        OverlapMode::ContainsRegion
    } else {
        OverlapMode::Overlap
    };

    let mut results = Vec::new();
    for (&seq_num, tree) in &index_data.chr_entries {
        let chr_regs = &buckets[seq_num as usize];
        if chr_regs.is_empty() { continue; }
        if index_data.common.verbose {
            eprintln!("[DEBUG] Querying chromosome {} with {} regions", seq_num, chr_regs.len());
        }
        for &(_, rstart, rend) in chr_regs {
            let hits: Vec<&Interval<u32>> = tree.query_interval(rstart, rend);
            if index_data.common.verbose {
                println!("Found {} hits for {} - {}.", hits.len(), rstart, rend);
            }
            results.extend(hits.into_iter().map(|iv| (iv.root_fid, iv.start, iv.end)));
        }
    }
    Ok(results)
}

/// Parse a single genomic region string (chr:start-end)
fn parse_region(region: &str, seqid_map: &FxHashMap<String, u32>, common: &CommonArgs) -> Result<(u32, u32, u32)> {
    let (seq, range) = region.split_once(':')
        .context("Invalid region format, expected 'chr:start-end'")?;
    let (s, e) = range.split_once('-')
        .context("Invalid range format, expected 'start-end'")?;
    let start = s.parse::<u32>()?;
    let end = e.parse::<u32>()?;
    let chr = seqid_map.get(seq)
        .with_context(|| format!("Sequence ID not found: {}", seq))?;
    if start >= end { anyhow::bail!("Region start must be less than end ({} >= {})", start, end); }
    if common.verbose {
        eprintln!("[DEBUG] Parsed region: chr={}, start={}, end={}", chr, start, end);
    }
    Ok((*chr, start, end))
}

/// Parse BED file using mmap + zero-copy field splitting
fn parse_bed_file(
    bed_path: &Path,
    seqid_map: &FxHashMap<String, u32>,
    common: &CommonArgs,
) -> Result<Vec<(u32, u32, u32)>> {
    let file = File::open(bed_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mut regions = Vec::new();
    for line in mmap.split(|&b| b == b'\n') {
        if line.is_empty() || line[0] == b'#' { continue; }
        let line_str = std::str::from_utf8(line)?;
        let mut parts = line_str.split_ascii_whitespace();
        if let (Some(seq), Some(s), Some(e)) = (parts.next(), parts.next(), parts.next()) {
            if let Some(&chr) = seqid_map.get(seq) {
                let start = parse::<u32>(s.as_bytes())?;
                let end = parse::<u32>(e.as_bytes())?;
                regions.push((chr, start, end));
            }
        }
    }
    if common.verbose {
        eprintln!("[DEBUG] Parsed {} regions from BED file", regions.len());
    }
    Ok(regions)
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
        let seqid_to_num = seqids.into_iter().enumerate().map(|(i, sq)| (sq, i as u32)).collect();
        let rit_path = append_suffix(index_prefix, ".rit");
        let rix_path = append_suffix(index_prefix, ".rix");
        let chr_entries = Self::load_region_index(&rit_path, &rix_path)?;
        Ok(Self { chr_entries, seqid_to_num, common: common.clone() })
    }

    fn load_region_index(
        rit_path: &Path,
        rix_path: &Path,
    ) -> Result<FxHashMap<u32, IntervalTree<u32>>>
    {
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

/// Main execution function
pub fn run(args: &IntersectArgs) -> Result<()> {
    let start_all = std::time::Instant::now();
    if args.common.verbose {
        eprintln!("[INFO] Loading regions from: {}", args.common.input.display());
    }

    let (_, seqid_map) = load_sqs(&args.common.input)?;
    let regions = if let Some(bed) = &args.bed {
        parse_bed_file(bed, &seqid_map, &args.common)?
    } else if let Some(r) = &args.region {
        vec![parse_region(r, &seqid_map, &args.common)?]
    } else {
        anyhow::bail!("No region specified");
    };

    let index_data = IndexData::load(&args.common.input, &args.common)?;
    let feats = query_features(
        &index_data,
        regions,
        args.contained,
        args.contains_region,
        args.invert,
    )?;

    let mut ids: Vec<u32> = feats.iter().map(|&( id, _, _)| id).collect();
    ids.sort_unstable(); ids.dedup();

    let gof = load_gof(&args.common.input)?;
    let gof_map: FxHashMap<_, _> = gof.into_iter()
        .map(|e| (e.feature_id, (e.start_offset, e.end_offset)))
        .collect();

    extract_gff_blocks(
        &args.common.input,
        &gof_map,
        &ids,
        &args.common.output,
        args.common.verbose,
    )?;

    if args.common.verbose {
        eprintln!("[TIME] Total: {:.2?}", start_all.elapsed());
    }
    Ok(())
}

/// Extract GFF blocks using memory-mapped file
fn extract_gff_blocks(
    gff_path: &Path,
    gof_map: &FxHashMap<u32, (u64, u64)>,
    feature_ids: &[u32],
    output: &Option<PathBuf>,
    verbose: bool,
) -> Result<()> {
    let file = File::open(gff_path)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let mut out: Box<dyn Write> = match output {
        Some(p) => Box::new(File::create(p)?),
        None    => Box::new(std::io::stdout()),
    };
    
    if verbose {
        eprintln!("[DEBUG] Extracting GFF blocks for {} features", feature_ids.len());
    }
    
    let mut extracted = 0;
    for &id in feature_ids {
        if let Some(&(s, e)) = gof_map.get(&id) {
            let (s, e) = (s as usize, e as usize);
            if s < e && e <= mmap.len() {
                out.write_all(&mmap[s..e])?;
                extracted += 1;
            } else if verbose {
                eprintln!("[WARNING] Invalid offsets for feature {}: {}..{}", id, s, e);
            }
        } else if verbose {
            eprintln!("[WARNING] Feature ID {} not found in GOF map", id);
        }
    }
    if verbose {
        eprintln!("[INFO] Successfully extracted {} features", extracted);
    }
    Ok(())
}
