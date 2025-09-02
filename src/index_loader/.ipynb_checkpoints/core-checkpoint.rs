use crate::append_suffix;
use anyhow::{Context, Result, bail};
use memmap2::Mmap;
use rustc_hash::FxHashMap;
use std::{
    //    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};



pub fn safe_mmap_readonly(path: &Path) -> Result<Mmap> {
    let file = File::open(path).with_context(|| format!("Failed to open file: {:?}", path))?;
    unsafe { Mmap::map(&file) }.with_context(|| format!("Failed to mmap file: {:?}", path))
}

pub fn load_sqs<P: AsRef<Path>>(path: P) -> Result<(Vec<String>, FxHashMap<String, u32>)> {
    let path = path.as_ref();
    let sqs_path = append_suffix(path, ".sqs");
    let file = File::open(&sqs_path)
        .with_context(|| format!("Failed to open SQS file: {:?}", &sqs_path))?;
    let reader = BufReader::new(file);

    let id_to_name: Vec<String> = reader.lines().collect::<Result<_, _>>()?;
    let name_to_id: FxHashMap<_, _> = id_to_name
        .iter()
        .enumerate()
        .map(|(id, name)| (name.clone(), id as u32))
        .collect();

    Ok((id_to_name, name_to_id))
}


pub fn load_atn(path: &Path) -> Result<(String, Vec<String>)> {
    let atn_path = append_suffix(path, ".atn");
    let mmap = safe_mmap_readonly(&atn_path)?;
    let data = &mmap[..];

    let mut values = Vec::new();
    let mut attr_name: Option<String> = None;

    // Helper to process a single line (without trailing '\n')
    let mut push_line = |bytes: &[u8]| -> Result<()> {
        if bytes.is_empty() {
            return Ok(());
        }
        let mut line = std::str::from_utf8(bytes)
            .context("ATN contains invalid UTF-8")?
            .trim();

        // Strip UTF-8 BOM if it appears at the start of the very first line
        if attr_name.is_none() && line.starts_with('\u{feff}') {
            line = &line['\u{feff}'.len_utf8()..];
        }

        if let Some(rest) = line.strip_prefix("#attribute=") {
            // Only a single header is allowed
            if attr_name.is_some() {
                bail!("Multiple #attribute= headers found in .atn file");
            }
            attr_name = Some(rest.to_string());
        } else if !line.is_empty() && !line.starts_with('#') {
            // Collect non-empty, non-comment value lines
            values.push(line.to_string());
        }
        Ok(())
    };

    // Scan by '\n'
    let mut start = 0usize;
    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            push_line(&data[start..i])?;
            start = i + 1;
        }
    }
    // Handle a trailing line without '\n'
    if start < data.len() {
        push_line(&data[start..])?;
    }

    let attr_name = attr_name
        .ok_or_else(|| anyhow::anyhow!("Missing #attribute=... header in .atn file"))?;

    Ok((attr_name, values))
}