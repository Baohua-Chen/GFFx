use anyhow::{Result, bail};
use clap::{ArgGroup, Parser};
use rayon::prelude::*;
use regex::Regex;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Write},
    path::PathBuf,
};

use crate::{
    CommonArgs, GofEntry, load_a2f, load_atn, load_fts, load_gof, load_prt, safe_mmap_readonly,
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

    #[arg(
        short = 'A',
        long,
        help = "Attribute list file (one per line)",
        group = "attr_input"
    )]
    attr_list: Option<PathBuf>,

    #[arg(
        short = 'a',
        long,
        help = "Single attribute value to search",
        group = "attr_input"
    )]
    attr: Option<String>,

    #[arg(short = 'r', long, help = "Enable regex mode for attribute matching")]
    regex: bool,
}

pub fn run(args: &SearchArgs) -> Result<()> {
    // 配置 Rayon 线程数（只设置一次全局线程池）
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.common.threads)
        .build_global()
        .ok(); // 已设置过时会返回 Err，这里忽略

    // 收集属性值：来自文件或单值
    let attr_values: Vec<String> = if let Some(file) = &args.attr_list {
        let reader = BufReader::new(File::open(file)?);
        reader
            .lines() // Iterator<Item = io::Result<String>>
            .map(|r| r.map(|s| s.trim().to_owned())) // trim
            .filter(|r| r.as_ref().map_or(true, |s| !s.is_empty()))
            .collect::<Result<Vec<_>, _>>()? // 遇到 Err 直接返回
    } else if let Some(val) = &args.attr {
        vec![val.clone()]
    } else {
        bail!("Either --attr-list (-A) or --attr (-a) must be provided.");
    };

    // 载入索引与原文件缓冲
    let gff_buf = safe_mmap_readonly(&args.common.input)?;
    let (atn_attr_name, atn_values) = load_atn(&args.common.input)?;
    let a2f = load_a2f(&args.common.input)?;
    let prt = load_prt(&args.common.input)?;
    let gof = load_gof(&args.common.input)?;
    let fts = load_fts(&args.common.input)?;

    // 索引一致性检查
    if a2f.len() != prt.len() {
        bail!(
            "Index error: a2f and prt length mismatch ({} vs {}). Possible feature registration inconsistency.",
            a2f.len(),
            prt.len()
        );
    }
    for (i, entry) in prt.iter().enumerate() {
        if (entry.parent as usize) >= prt.len() {
            bail!(
                "Index error: feature {} refers to non-existent parent {}",
                i,
                entry.parent
            );
        }
    }

    // 快速按 feature_id 查 GOF
    let gof_map: HashMap<u32, GofEntry> = gof.into_iter().map(|e| (e.feature_id, e)).collect();

    // 将“待匹配的属性值”映射到其 attr_id 列表
    let mut attr_to_ids: HashMap<String, Vec<u32>> = HashMap::new();

    if args.regex {
        // 正确把 Vec<String> -> Vec<Regex>
        let patterns: Vec<Regex> = attr_values
            .iter()
            .map(String::as_str) // &String -> &str
            .map(Regex::new)
            .collect::<std::result::Result<Vec<_>, _>>()?;

        for (i, val) in atn_values.iter().enumerate() {
            if patterns.iter().any(|re| re.is_match(val)) {
                attr_to_ids.entry(val.clone()).or_default().push(i as u32);
            }
        }
    } else {
        // O(1) 集合查找
        let wanted: HashSet<&str> = attr_values.iter().map(String::as_str).collect();
        for (i, val) in atn_values.iter().enumerate() {
            if wanted.contains(val.as_str()) {
                attr_to_ids.entry(val.clone()).or_default().push(i as u32);
            }
        }
    }

    if attr_to_ids.is_empty() {
        bail!(
            "None of the attributes matched. Index was built with attribute: {}",
            atn_attr_name
        );
    }

    // 并行构造输出块
    let results: Vec<String> = attr_to_ids
        .par_iter()
        .map(|(attr, attr_ids)| {
            let mut roots_seen = HashSet::new();
            let mut lines = Vec::new();

            for &aid in attr_ids {
                // 遍历 a2f 上所有具有该 attr_id 的要素
                for a2f_entry in a2f.iter().filter(|e| e.attr_id == aid) {
                    let matched_fid = a2f_entry.feature_id;

                    // 找祖先根（parent == self）
                    let mut fid = matched_fid;
                    while prt[fid as usize].parent != fid {
                        fid = prt[fid as usize].parent;
                    }

                    // 重复的根模型跳过
                    if !roots_seen.insert(fid) {
                        continue;
                    }

                    // 取该根模型对应 GOF 区块，没有就跳过
                    let Some(entry) = gof_map.get(&fid) else {
                        continue;
                    };

                    let slice = &gff_buf[entry.start_offset as usize..entry.end_offset as usize];
                    let mut block = format!(
                        "# match attribute: {} (via feature_id={} in model={})\n",
                        attr,
                        matched_fid,
                        fts.get(fid as usize).unwrap_or(&"<unknown>".to_string())
                    );
                    let text = std::str::from_utf8(slice)
                        .map_err(|_| anyhow::anyhow!("Invalid UTF-8 in extracted block"))?;
                    block.push_str(text);
                    lines.push(block);
                }
            }

            Ok::<_, anyhow::Error>(lines.join(""))
        })
        .collect::<Result<Vec<_>>>()?;

    // 输出
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
