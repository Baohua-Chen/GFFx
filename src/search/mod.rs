use anyhow::{Result, bail};
use clap::{Parser, ArgGroup};
use rayon::prelude::*;
use regex::Regex;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufReader, BufRead, Write},
    path::PathBuf,
};

use crate::{
    load_atn, load_a2f, load_prt, load_gof, safe_mmap_readonly, load_fts,
    CommonArgs, GofEntry,
};

#[derive(Parser, Debug)]
#[command(
    about = "Search features by attribute values and extract their full models.",
    group = ArgGroup::new("attr_input")
        .required(true)
        .args(["attr_list", "attr"])
)]
pub struct SearchArgs {
    /// Common input/output/thread arguments
    #[clap(flatten)]
    pub common: CommonArgs,

    #[arg(short = 'A', long, help = "Attribute list file (one per line)", group = "attr_input")]
    attr_list: Option<PathBuf>,

    #[arg(short = 'a', long, help = "Single attribute value to search", group = "attr_input")]
    attr: Option<String>,

    #[arg(short = 'r', long, help = "Enable regex mode for attribute matching")]
    regex: bool,
}

pub fn run(args: &SearchArgs) -> Result<()> {
    rayon::ThreadPoolBuilder::new().num_threads(args.common.threads).build_global()?;

    let attr_values: Vec<String> = if let Some(file) = &args.attr_list {
        BufReader::new(File::open(file)?).lines().collect::<Result<_, _>>()?
    } else if let Some(val) = &args.attr {
        vec![val.clone()]
    } else {
        bail!("Either --attr-list (-A) or --attr (-a) must be provided.");
    };

    let gff_buf = safe_mmap_readonly(&args.common.input)?;
    let (atn_attr_name, atn_values) = load_atn(&args.common.input)?;
    let a2f = load_a2f(&args.common.input)?;
    let prt = load_prt(&args.common.input)?;
    let gof = load_gof(&args.common.input)?;
    let fts = load_fts(&args.common.input)?;

    if a2f.len() != prt.len() {
        bail!("Index error: a2f and prt length mismatch ({} vs {}). Possible feature registration inconsistency.", a2f.len(), prt.len());
    }

    for (i, entry) in prt.iter().enumerate() {
        if (entry.parent as usize) >= prt.len() {
            bail!("Index error: feature {} refers to non-existent parent {}", i, entry.parent);
        }
    }

    let gof_map: HashMap<u32, GofEntry> = gof.into_iter().map(|e| (e.feature_id, e)).collect();

    let mut attr_to_ids: HashMap<String, Vec<u32>> = HashMap::new();

    if args.regex {
        let patterns: Vec<Regex> = attr_values.iter()
            .map(|pat| Regex::new(pat))
            .collect::<std::result::Result<_, _>>()?;

        for (i, val) in atn_values.iter().enumerate() {
            if patterns.iter().any(|re| re.is_match(val)) {
                attr_to_ids.entry(val.clone()).or_default().push(i as u32);
            }
        }
    } else {
        for (i, val) in atn_values.iter().enumerate() {
            if attr_values.contains(val) {
                attr_to_ids.entry(val.clone()).or_default().push(i as u32);
            }
        }
    }

    if attr_to_ids.is_empty() {
        bail!("{}", format!("None of the attributes matched. Index was built with attribute: {}", atn_attr_name));
    }

    let results: Vec<String> = attr_to_ids.par_iter()
        .map(|(attr, attr_ids)| {
            let mut roots_seen = HashSet::new();
            let mut lines = Vec::new();

            for &aid in attr_ids {
                for a2f_entry in a2f.iter().filter(|e| e.attr_id == aid) {
                    let matched_fid = a2f_entry.feature_id;
                    let mut fid = matched_fid;
                    while prt[fid as usize].parent != fid {
                        fid = prt[fid as usize].parent;
                    }
                    if roots_seen.insert(fid) {
                        if let Some(entry) = gof_map.get(&fid) {
                            let slice = &gff_buf[entry.start_offset as usize..entry.end_offset as usize];
                            let mut block = format!(
                                "# match attribute: {} (via feature_id={} in model={})\n",
                                attr, matched_fid, fts.get(fid as usize).unwrap_or(&"<unknown>".to_string())
                            );
                            let text = std::str::from_utf8(slice)
                                .map_err(|_| anyhow::anyhow!("Invalid UTF-8 in extracted block"))?;
                            block.push_str(text);
                            lines.push(block);
                        }
                    }
                }
            }

            Ok::<_, anyhow::Error>(lines.join(""))
        })
        .collect::<Result<Vec<_>>>()?;

    let out: Box<dyn Write> = match &args.common.output {
        Some(path) => Box::new(File::create(path)?),
        None => Box::new(std::io::stdout()),
    };
    let mut out = std::io::BufWriter::new(out);

    for block in results {
        out.write_all(block.as_bytes())?;
    }

    Ok(())
}

