use crate::append_suffix;
use crate::{Interval, IntervalTree, save_multiple_trees, write_offsets_to_file};
use anyhow::anyhow;
use anyhow::{Result, bail};
use byteorder::{LittleEndian, WriteBytesExt};
use indexmap::IndexMap;
use memchr::memchr;
use memmap2::Mmap;
use regex::{Regex, escape};
use std::{collections::HashMap, fs::File, io::Write, path::PathBuf};

// Writes text lines to a file
pub fn write_lines(path: PathBuf, lines: &[String]) -> Result<()> {
    let mut file = File::create(path)?;
    for line in lines {
        writeln!(file, "{}", line)?;
    }
    Ok(())
}

// Writes u32 values in binary little-endian format
pub fn write_binary_u32(path: PathBuf, values: &[u32]) -> Result<()> {
    let mut file = File::create(path)?;
    for &v in values {
        file.write_u32::<LittleEndian>(v)?;
    }
    Ok(())
}

// Writes GFF offset records (gof)
pub fn write_gof(file: &mut File, id: u32, start: u64, end: u64) -> Result<()> {
    file.write_u32::<LittleEndian>(id)?;
    file.write_u32::<LittleEndian>(0)?; // Reserved field
    file.write_u64::<LittleEndian>(start)?;
    file.write_u64::<LittleEndian>(end)?;
    Ok(())
}

/// Builds various index files for a GFF: .fts, .prt, .a2f, .atn, .sqs, .gof, .rit, .rix
pub fn build_index(gff: &PathBuf, attr_key: &str, verbose: bool) -> Result<()> {
    // Compile regex patterns
    let id_re = Regex::new(r"ID=([^;\s]+)")?;
    let parent_re = Regex::new(r"Parent=([^;\s]+)")?;
    let attr_re = Regex::new(&format!(r"{}=([^;\s]+)", escape(attr_key)))?;

    if verbose {
        eprintln!("Building index for {} ...", gff.display());
    }

    // Memory-map input file
    let file = File::open(gff)?;
    let mmap = unsafe { Mmap::map(&file)? };
    let data = &mmap[..];

    // First pass: parse raw features
    struct RawFeature {
        seqid: String,
        start: u32,
        end: u32,
        line_offset: u64,
        id: String,
        parent: Option<String>,
        attr: Option<String>,
    }
    let mut raw_features = Vec::new();
    let mut offset = 0;

    while offset < data.len() {
        let nl_pos = memchr(b'\n', &data[offset..])
            .map(|pos| pos + offset)
            .unwrap_or(data.len());
        let line_bytes = &data[offset..nl_pos];
        let line_offset = offset as u64;
        offset = nl_pos + 1;

        if line_bytes.is_empty() || line_bytes[0] == b'#' {
            continue;
        }
        let line = std::str::from_utf8(line_bytes)?.trim();
        if line.is_empty() {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            bail!("Invalid GFF line (expected 9 columns): {}", line);
        }

        let seqid = fields[0].to_string();
        let start = fields[3].parse::<u32>()?;
        let end = fields[4].parse::<u32>()?;

        // Extract ID
        let id = id_re
            .captures(line)
            .ok_or_else(|| anyhow!("Missing ID in feature: {}", line))?[1]
            .to_string();
        // Extract raw Parent (may refer to unseen ID)
        let parent = parent_re.captures(line).map(|cap| cap[1].to_string());
        // Extract attribute value
        let attr = attr_re.captures(line).map(|cap| cap[1].to_string());

        raw_features.push(RawFeature {
            seqid,
            start,
            end,
            line_offset,
            id,
            parent,
            attr,
        });
    }

    // Build feature_map: string ID -> numeric ID
    let mut feature_map: HashMap<String, u32> = HashMap::new();
    for (i, rf) in raw_features.iter().enumerate() {
        feature_map.insert(rf.id.clone(), i as u32);
    }

    // Open output files
    let mut fts_file = File::create(append_suffix(gff, ".fts"))?;
    let mut prt_entries = Vec::with_capacity(raw_features.len());
    let mut a2f_entries = Vec::with_capacity(raw_features.len());
    let mut atn_entries = Vec::new();
    let mut attr_value_to_id: HashMap<String, u32> = HashMap::new();
    let mut gof_file = File::create(append_suffix(gff, ".gof"))?;
    let mut seqid_intervals: IndexMap<String, Vec<(u32, u32, u32)>> = IndexMap::new();
    let mut current_root: Option<(u32, u64)> = None;

    // Write .fts and build .prt, .a2f, .gof, and seqid intervals
    for rf in &raw_features {
        let fid = feature_map[&rf.id];
        writeln!(fts_file, "{}", rf.id)?;
        // Resolve parent (fallback to self if missing)
        let parent_id = rf
            .parent
            .as_ref()
            .and_then(|p| feature_map.get(p).cloned())
            .unwrap_or(fid);
        prt_entries.push(parent_id);
        // Record roots for GOF and intervals
        if parent_id == fid {
            seqid_intervals
                .entry(rf.seqid.clone())
                .or_default()
                .push((rf.start, rf.end, fid));
            if let Some((old_id, old_off)) = current_root.take() {
                write_gof(&mut gof_file, old_id, old_off, rf.line_offset)?;
            }
            current_root = Some((fid, rf.line_offset));
        }
        // Attribute mapping
        if let Some(val) = &rf.attr {
            let aid = *attr_value_to_id.entry(val.clone()).or_insert_with(|| {
                let a = atn_entries.len() as u32;
                atn_entries.push(val.clone());
                a
            });
            a2f_entries.push(aid);
        } else {
            a2f_entries.push(u32::MAX);
        }
    }
    // Write final GOF record
    if let Some((last_id, last_off)) = current_root {
        write_gof(&mut gof_file, last_id, last_off, data.len() as u64)?;
    }

    // Build interval trees per seqid
    let mut seqid_to_num: IndexMap<String, u32> = IndexMap::new();
    let mut trees = Vec::with_capacity(seqid_intervals.len());
    for (seqid, ivs) in seqid_intervals.iter() {
        let num = seqid_to_num.len() as u32;
        seqid_to_num.insert(seqid.clone(), num);
        let iv_structs: Vec<Interval<_>> = ivs
            .iter()
            .map(|&(start, end, fid)| Interval {
                start,
                end,
                root_fid: fid,
            })
            .collect();
        trees.push(IntervalTree::new(iv_structs));
    }

    // Write .rit and .rix
    let rit = append_suffix(gff, ".rit");
    let rix = append_suffix(gff, ".rix");
    let offsets = save_multiple_trees(rit.to_str().unwrap(), &trees)?;
    write_offsets_to_file(&offsets, rix.to_str().unwrap())?;

    // Write .sqs (sequence list)
    let seqids: Vec<String> = seqid_to_num.keys().cloned().collect();
    write_lines(append_suffix(gff, ".sqs"), &seqids)?;

    // Write .atn, .a2f, .prt
    let mut atn_out = Vec::with_capacity(atn_entries.len() + 1);
    atn_out.push(format!("#attribute={}", attr_key));
    atn_out.extend(atn_entries.clone());
    write_lines(append_suffix(gff, ".atn"), &atn_out)?;
    write_binary_u32(append_suffix(gff, ".a2f"), &a2f_entries)?;
    write_binary_u32(append_suffix(gff, ".prt"), &prt_entries)?;

    if verbose {
        eprintln!("Index built successfully for {}", gff.display());
    }
    Ok(())
}
