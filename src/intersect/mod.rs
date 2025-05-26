use anyhow::{Result, Context, bail};
use clap::{Parser, ArgGroup};
use memmap2::MmapOptions;
use rayon::{prelude::*, ThreadPoolBuilder};

use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Write, stdout},
    path::PathBuf,
};

use crate::{
    CommonArgs, load_sqs, load_gbi, load_gof, GbiEntry,
    write_gff_header, write_block
};

#[derive(Parser, Debug)]
#[command(
    about = "Extract models by a region or regions from a BED file",
    long_about = "This tool extracts features and their parent models intersect with a specified region or regions from a BED file"
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

    #[arg(short = 'r', long, group = "regions")]
    pub region: Option<String>,

    #[arg(short = 'b', long, group = "regions")]
    pub bed: Option<PathBuf>,

    #[arg(short = 'c', long, group = "mode")]
    pub contained: bool,

    #[arg(short = 'C', long, group = "mode")]
    pub contains_region: bool,

    #[arg(short = 'O', long, group = "mode")]
    pub overlap: bool,

    #[arg(short = 'v', long, default_value_t = false)]
    pub invert: bool,
}

#[derive(Debug)]
pub struct GbiBinIndex {
    pub bin_map: HashMap<(u32, u32), Vec<(u32, u32, u32)>>, // (feature_id, start, end)
}

impl GbiBinIndex {
    pub fn from_entries(entries: &[GbiEntry]) -> Self {
        let mut bin_map: HashMap<(u32, u32), Vec<(u32, u32, u32)>> = HashMap::new();
        for e in entries {
            bin_map.entry((e.seqid_id, e.bin))
                .or_default()
                .push((e.feature_id, e.start, e.end));
        }
        GbiBinIndex { bin_map }
    }

    pub fn query(&self, seqid_id: u32, start: u32, end: u32) -> Vec<(u32, u32, u32)> {
        let mut seen = HashSet::new();
        let mut results = Vec::new();
        for bin in reg2bins(start, end) {
            if let Some(entries) = self.bin_map.get(&(seqid_id, bin)) {
                for &(fid, f_start, f_end) in entries {
                    if seen.insert(fid) {
                        results.push((fid, f_start, f_end));
                    }
                }
            }
        }
        results
    }
}

fn reg2bins(start: u32, end: u32) -> Vec<u32> {
    let end = end - 1;
    let mut bins = vec![];
    for (shift, offset) in &[(17, 0), (20, 1), (23, 9), (26, 73), (29, 585)] {
        for i in (start >> shift)..=(end >> shift) {
            bins.push(offset + i);
        }
    }
    bins
}

fn parse_region(r: &str) -> Result<(String, u32, u32)> {
    let (seq, range) = r.split_once(':')
        .context("Region must be in 'seqid:start-end' format")?;
    let mut parts = range.splitn(2, '-');
    let start = parts.next().unwrap_or("1").parse().unwrap_or(1);
    let end = parts.next().unwrap_or("").parse().unwrap_or(u32::MAX);
    Ok((seq.to_string(), start, end))
}

fn load_bed(path: &PathBuf) -> Result<Vec<(String, String)>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut regs = Vec::new();
    for line in reader.lines() {
        let l = line?;
        if l.starts_with('#') || l.trim().is_empty() { continue; }
        let cols: Vec<&str> = l.split_whitespace().collect();
        if cols.len() < 3 { continue; }
        let start = cols[1].parse::<u32>().unwrap_or(0) + 1;
        let end = cols[2].parse::<u32>().unwrap_or(start);
        let spec = format!("{}:{}-{}", cols[0], start, end);
        regs.push((spec, l));
    }
    Ok(regs)
}

pub fn run(args: &IntersectArgs) -> Result<()> {
    let _ = ThreadPoolBuilder::new().num_threads(args.common.threads).build_global();

    if [args.contained, args.contains_region, args.overlap].iter().filter(|&&x| x).count() > 1 {
        bail!("Only one of --contained, --contains-region, or --overlap can be set.");
    }

    let gff_path = &args.common.input;
    let file = File::open(gff_path)?;
    let mmap = unsafe { MmapOptions::new().map(&file)? };
    let gff_buf = &mmap[..];

    let seqs = load_sqs(gff_path)?;
    let seq_map: HashMap<String, u32> = seqs.iter().enumerate()
        .map(|(i, s)| (s.clone(), i as u32)).collect();
    let gof_map: HashMap<u32, (usize, usize)> = load_gof(gff_path)?
        .into_iter().map(|e| (e.feature_id, (e.start_offset as usize, e.end_offset as usize))).collect();
    let gbi_index = GbiBinIndex::from_entries(&load_gbi(gff_path)?);

    let mut regions = Vec::new();
    if let Some(r) = &args.region { regions.push((r.clone(), format!("Region: {}", r))); }
    if let Some(bf) = &args.bed { regions.extend(load_bed(bf)?); }
    if regions.is_empty() { bail!("No region or BED specified"); }

    let results: Vec<_> = regions.par_iter().filter_map(|(spec, label)| {
        let (seq, rs, re) = parse_region(spec).ok()?;
        let &seqid = seq_map.get(&seq)?;
        let entries = gbi_index.query(seqid, rs, re);
        let offsets = entries.into_iter().filter(|&(_, s, e)| {
            let sel = if args.contains_region {
                s <= rs && e >= re
            } else if args.contained {
                s >= rs && e <= re
            } else {
                s <= re && e >= rs
            };
            if args.invert { !sel } else { sel }
        }).filter_map(|(fid, _, _)| gof_map.get(&fid).copied()).collect::<Vec<_>>();
        Some((label.clone(), offsets))
    }).collect();

    let mut writer: Box<dyn Write> = match &args.common.out {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(stdout()),
    };

    write_gff_header(&mut writer, gff_buf)?;
    for (label, blocks) in results {
        writeln!(writer, "# {}", label)?;
        for (so, eo) in blocks {
            write_block(&mut writer, &gff_buf[so..eo], args.common.types.as_deref())?;
        }
    }
    Ok(())
}
