use anyhow::{Result, bail};
use byteorder::{LittleEndian, WriteBytesExt};
use regex::Regex;
use crate::append_suffix;
use std::{
    collections::{HashMap, BTreeMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Seek, SeekFrom, Write},
    path::{Path, PathBuf},
};

fn write_lines<P: AsRef<Path>>(path: P, lines: &[String]) -> Result<()> {
    let mut file = File::create(path)?;
    for line in lines {
        writeln!(file, "{}", line)?;
    }
    Ok(())
}

fn write_binary_u32<P: AsRef<Path>>(path: P, data: &[u32]) -> Result<()> {
    let mut file = File::create(path)?;
    for v in data {
        file.write_u32::<LittleEndian>(*v)?;
    }
    Ok(())
}

pub fn build_index(gff: &PathBuf, attr_key: &str) -> Result<()> {
    let id_re = Regex::new(r"ID=([^;\s]+)")?;
    let parent_re = Regex::new(r"Parent=([^;\s]+)")?;
    let attr_re = Regex::new(&format!(r"{}=([^;\s]+)", regex::escape(attr_key)))?;

    let mut reader = BufReader::new(File::open(&gff)?);
    let mut prev_offset: u64 = 0;

    let mut feature_counter: u32 = 0;
    let mut feature_map: HashMap<String, u32> = HashMap::new(); // ID -> feature index
    let mut a2f_entries: Vec<u32> = Vec::new();                 // index -> attr_id
    let mut prt_entries: Vec<u32> = Vec::new();                 // index -> parent index

    let mut attr_value_to_id: HashMap<String, u32> = HashMap::new(); // attr val -> attr id
    let mut atn_entries: Vec<String> = Vec::new();                   // ordered attr val

    let mut seqid_set: HashSet<String> = HashSet::new();
    let mut seqid_list: Vec<String> = Vec::new();

    let mut root_features: Vec<(u32, String, u32, u32)> = Vec::new(); // for .gbi
    let mut current_root: Option<(u32, u64)> = None; // (id, start offset)

    let mut gof_file = File::create(append_suffix(&gff, ".gof"))?;
    let mut fts_file = File::create(append_suffix(&gff, ".fts"))?;

    loop {
        let offset = reader.seek(SeekFrom::Current(0))?;
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            if let Some((id, start)) = current_root.take() {
                write_gof(&mut gof_file, id, start, prev_offset)?;
            }
            break;
        }

        prev_offset = offset + line.len() as u64;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 9 {
            bail!("GFF line does not contain 9 columns: {}", line);
        }

        let (seqid, start, end) = (fields[0], fields[3].parse::<u32>()?, fields[4].parse::<u32>()?);
        if seqid_set.insert(seqid.to_string()) {
            seqid_list.push(seqid.to_string());
        }

        // === ID ===
        let id = if let Some(cap) = id_re.captures(line) {
            cap[1].to_string()
        } else {
            bail!("Feature missing ID: {}", line);
        };

        let (feature_id, is_new) = if let Some(&fid) = feature_map.get(&id) {
            (fid, false)
        } else {
            let fid = feature_counter;
            feature_map.insert(id.clone(), fid);
            feature_counter += 1;
            writeln!(fts_file, "{}", id)?;
            (fid, true)
        };

        // === Parent ===
        let parent_id = if let Some(cap) = parent_re.captures(line) {
            let parent_str = &cap[1];
            match feature_map.get(parent_str) {
                Some(&pid) => pid,
                None => bail!("Parent '{}' not found before child '{}'", parent_str, id),
            }
        } else {
            feature_id // self parent
        };

        if is_new {
            prt_entries.push(parent_id);

            if parent_id == feature_id {
                if let Some((old_id, old_start)) = current_root.take() {
                    write_gof(&mut gof_file, old_id, old_start, offset)?;
                }
                current_root = Some((feature_id, offset));
                root_features.push((feature_id, seqid.to_string(), start, end));
            }
        } else {
            // check parent consistency
            let old_pid = prt_entries[feature_id as usize];
            if old_pid != parent_id {
                bail!("Conflicting Parent for ID={}: {} vs {}", id, old_pid, parent_id);
            }
        }

        // === Attribute ===
        let attr_val = if let Some(cap) = attr_re.captures(line) {
            Some(cap[1].to_string())
        } else {
            None
        };

        if is_new {
            let attr_id = if let Some(v) = attr_val {
                *attr_value_to_id.entry(v.clone()).or_insert_with(|| {
                    let aid = atn_entries.len() as u32;
                    atn_entries.push(v);
                    aid
                })
            } else {
                u32::MAX
            };
            a2f_entries.push(attr_id);
        } else if let Some(v) = attr_val {
            let attr_id = *attr_value_to_id.entry(v.clone()).or_insert_with(|| {
                let aid = atn_entries.len() as u32;
                atn_entries.push(v);
                aid
            });
            let old_id = a2f_entries[feature_id as usize];
            if old_id != attr_id {
                bail!("Conflicting {} for ID={}: {} vs {}", attr_key, id, old_id, attr_id);
            }
        }
    }

    if atn_entries.is_empty() {
        eprintln!(
            "Warning: attribute '{}' not found in any features. \
            Search by attribute will not be available.",
            attr_key
        );
    }

    write_lines(append_suffix(gff, ".sqs"), &seqid_list)?;
    write_lines(append_suffix(gff, ".atn"), &{
        let mut lines = vec![format!("#attribute={}", attr_key)];
        lines.extend(atn_entries.iter().cloned());
        lines
    })?;
    write_binary_u32(append_suffix(gff, ".a2f"), &a2f_entries)?;
    write_binary_u32(append_suffix(gff, ".prt"), &prt_entries)?;
    write_gbi(&append_suffix(gff, ".gbi"), &root_features)?;
    Ok(())
}

fn write_gof(file: &mut File, id: u32, start: u64, end: u64) -> Result<()> {
    file.write_u32::<LittleEndian>(id)?;
    file.write_u32::<LittleEndian>(0)?;
    file.write_u64::<LittleEndian>(start)?;
    file.write_u64::<LittleEndian>(end)?;
    Ok(())
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

fn write_gbi(path: &Path, features: &[(u32, String, u32, u32)]) -> Result<()> {
    let mut file = File::create(path)?;
    let mut seqid_map: HashMap<String, u32> = HashMap::new();
    let mut seqid_list: Vec<String> = Vec::new();
    let mut bin_map: BTreeMap<(u32, u32), Vec<(u32, u32, u32)>> = BTreeMap::new();

    for &(id, ref seqid, start, end) in features {
        let seqid_id = *seqid_map.entry(seqid.clone()).or_insert_with(|| {
            let id = seqid_list.len() as u32;
            seqid_list.push(seqid.clone());
            id
        });
        for bin in reg2bins(start, end) {
            bin_map.entry((seqid_id, bin)).or_default().push((id, start, end));
        }
    }

    file.write_u32::<LittleEndian>(seqid_list.len() as u32)?;
    let offset_pos = file.stream_position()?;
    for _ in &seqid_list {
        file.write_u64::<LittleEndian>(0)?;
    }

    let mut offsets = vec![0u64; seqid_list.len()];
    for (i, _) in seqid_list.iter().enumerate() {
        let pos = file.stream_position()?;
        offsets[i] = pos;
        let bins: Vec<_> = bin_map.iter().filter(|((sid, _), _)| *sid == i as u32).collect();

        file.write_u32::<LittleEndian>(i as u32)?;
        file.write_u32::<LittleEndian>(bins.len() as u32)?;
        for ((_, bin), entries) in bins {
            file.write_u32::<LittleEndian>(*bin)?;
            file.write_u32::<LittleEndian>(entries.len() as u32)?;
            for &(id, start, end) in entries {
                file.write_u32::<LittleEndian>(id)?;
                file.write_u32::<LittleEndian>(start)?;
                file.write_u32::<LittleEndian>(end)?;
            }
        }
    }

    file.seek(SeekFrom::Start(offset_pos))?;
    for offset in offsets {
        file.write_u64::<LittleEndian>(offset)?;
    }
    Ok(())
}

pub fn check_index_files_exist(gff: &PathBuf, rebuild: bool, attr_key: &str) -> Result<bool> {
    let expected_suffixes = [".gof", ".fts", ".prt", ".gbi", ".sqs", ".atn", ".a2f"];
    let mut missing = Vec::new();

    for ext in &expected_suffixes {
        let path = append_suffix(gff, ext);
        if !path.exists() {
            missing.push(ext.to_string());
        }
    }

    if !missing.is_empty() {
        let joined = missing.join(", ");
        eprintln!("Some index files are missing. Rebuilding is required: {}", joined);

        if rebuild {
            build_index(gff, attr_key)?;
            return Ok(true);
        }

        return Ok(false);
    }

    Ok(true)
}
