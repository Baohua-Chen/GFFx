use anyhow::{Result, Context, bail};
use clap::Parser;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufReader, BufRead},
    path::PathBuf,
    sync::{Arc, Mutex},
};
use rayon::prelude::*;
use crate::{
    CommonArgs, load_fts, load_prt, load_gof, 
    write_gff_output,
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

    #[arg(short = 'f', long, group = "feature")]
    pub feature_id: Option<String>,

    #[arg(short = 'F', long, group = "feature")]
    pub feature_file: Option<PathBuf>,

    #[arg(short = 'd', long)]
    pub descendants_only: bool,
}

/// Recursively walk up parent map to find root feature ID, with loop detection
fn resolve_root(start: u32, prt_map: &HashMap<u32, u32>) -> Result<u32> {
    let mut current = start;
    let mut visited = HashSet::new();
    while visited.insert(current) {
        match prt_map.get(&current) {
            Some(&parent) if parent == current => return Ok(current),
            Some(&parent) => current = parent,
            None => bail!("Parent not found for feature ID {}", current),
        }
    }
    bail!("Circular parent relationship detected for feature ID {}", start)
}

fn extract_features(
    feature_ids: &HashSet<String>,
    fts_map: &HashMap<&str, u32>,
    prt_map: &HashMap<u32, u32>,
    gof_map: &HashMap<u32, (u64, u64)>,
    descendants_only: bool,
    threads: usize,
) -> Vec<(u64, u64)> {
    let seen_roots: Arc<Mutex<HashSet<u32>>> = Arc::new(Mutex::new(HashSet::new()));
    let process_feature = |fid: &String| -> Option<(u64, u64)> {
        let numeric_id = *fts_map.get(fid.as_str())?;
        let root_id = resolve_root(numeric_id, prt_map).ok()?;

        {
            let mut seen = seen_roots.lock().unwrap();
            if !seen.insert(root_id) {
                return None;
            }
        }

        if descendants_only {
            eprintln!("Descendants only mode will implemented in the next version.");
            return None;
        }

        gof_map.get(&root_id).cloned()
    };

/*
            match process_descendants(fid, fts_map, prt_map, gof_map, gff_path) {
                Ok(slice) => {
                    let mut file = tempfile::tempfile().ok()?;
                    use std::io::Write;
                    file.write_all(&slice).ok()?;
                    Some((0, slice.len() as u64))
                },
                Err(err) => {
                    eprintln!("Warning: {} (feature '{}')", err, fid);
                    None
                }
            }
        } else {
            gof_map.get(&root_id).copied()
        }
    };
*/
    if threads > 1 && feature_ids.len() > 2 {
        feature_ids.par_iter().filter_map(process_feature).collect()
    } else {
        feature_ids.iter().filter_map(process_feature).collect()
    }
}

pub fn run(args: &ExtractArgs) -> Result<()> {
    let gff_path = &args.common.input;
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.common.threads)
        .build_global()
        .expect("Failed to set thread pool size");

    let fts = load_fts(gff_path)?;
    let prt = load_prt(gff_path)?;
    let gof = load_gof(gff_path)?;

    let fts_map: HashMap<&str, u32> = fts.iter().enumerate().map(|(i, s)| (s.as_str(), i as u32)).collect();
    let prt_map: HashMap<u32, u32> = prt.iter().map(|e| (e.child, e.parent)).collect();
    let gof_map: HashMap<u32, (u64, u64)> = gof.iter().map(|e| (e.feature_id, (e.start_offset, e.end_offset))).collect();

    let feature_ids: HashSet<String> = if let Some(ref file_path) = args.feature_file {
        let reader = BufReader::new(File::open(file_path)
            .with_context(|| format!("Cannot open feature list: {:?}", file_path))?);
        reader
            .lines()
            .filter_map(Result::ok)
            .filter(|s| !s.trim().is_empty())
            .collect()
    } else if let Some(ref single_id) = args.feature_id {
        [single_id.clone()].iter().cloned().collect()
    } else {
        bail!("Either --feature-id (-f) or --feature-file (-F) must be specified");
    };

    // Track seen root feature IDs (to avoid duplication)
    let block_offsets = extract_features(
        &feature_ids,
        &fts_map,
        &prt_map,
        &gof_map,
        args.descendants_only,
        args.common.threads,
    );

    write_gff_output(
        gff_path,
        &block_offsets,
        &args.common.output,
        args.common.types.as_deref(),
        args.common.verbose,
    )?;

    Ok(())
}

/*
fn process_one<'a>(
    fid: &str,
    fts_map: &HashMap<&str, u32>,
    prt_map: &HashMap<u32, u32>,
    gof_map: &HashMap<u32, (u64, u64)>,
    gff_buf: &'a [u8],
    descendants_only: bool,
) -> Result<Vec<u8>> {
    let numeric_id = *fts_map
        .get(fid)
        .with_context(|| format!("Feature not found in .fts"))?;

    let root_id = resolve_root(numeric_id, prt_map)?;
    let (start, end) = *gof_map
        .get(&root_id)
        .with_context(|| format!("Root feature {} not found in .gof", root_id))?;
    let slice = &gff_buf[start as usize..end as usize];

    if !descendants_only {
        return Ok(slice.to_vec());
    }

    let text = std::str::from_utf8(slice)?.trim_end();
    let lines: Vec<&str> = text.lines().collect();
    let mut id_to_children: HashMap<String, Vec<String>> = HashMap::new();
    let mut id_to_line: HashMap<String, &str> = HashMap::new();

    for line in &lines {
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        if let Some(attrs) = line.split('\t').nth(8) {
            let mut id: Option<String> = None;
            let mut parents: Vec<String> = vec![];
            for kv in attrs.split(';') {
                let kv = kv.trim();
                if kv.starts_with("ID=") {
                    id = Some(kv[3..].to_string());
                } else if kv.starts_with("Parent=") {
                    parents = kv[7..].split(',').map(|s| s.to_string()).collect();
                }
            }
            if let Some(id) = id {
                id_to_line.insert(id.clone(), line);
                for p in parents {
                    id_to_children.entry(p).or_default().push(id.clone());
                }
            }
        }
    }

    let mut keep_ids = HashSet::new();
    let mut stack = vec![fid.to_string()];
    while let Some(id) = stack.pop() {
        if keep_ids.insert(id.clone()) {
            if let Some(children) = id_to_children.get(&id) {
                stack.extend_from_slice(children);
            }
        }
    }

    let mut result = Vec::new();
    for id in &keep_ids {
        if let Some(line) = id_to_line.get(id) {
            result.push(*line);
        }
    }
    let result_str = result.join("\n") + "\n";
    Ok(result_str.into_bytes())
}
*/

