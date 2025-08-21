use crate::append_suffix;
use anyhow::{Context, Result, bail};
use byteorder::{ByteOrder, LittleEndian, ReadBytesExt};
use memmap2::Mmap;
use rustc_hash::FxHashMap;
use std::{
    //    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, ErrorKind},
    path::Path,
};

#[derive(Debug)]
pub struct GofEntry {
    pub feature_id: u32,
    pub start_offset: u64,
    pub end_offset: u64,
}

#[derive(Debug, Clone, Copy)]
pub struct PrtEntry {
    pub child: u32,
    pub parent: u32,
}

#[derive(Debug)]
pub struct GbiEntry {
    pub seqid_id: u32,
    pub bin: u32,
    pub feature_id: u32,
    pub start: u32,
    pub end: u32,
}

#[derive(Debug)]
pub struct A2fEntry {
    pub feature_id: u32,
    pub attr_id: u32,
}

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

pub fn load_gof<P: AsRef<Path>>(base: P) -> Result<Vec<GofEntry>> {
    let path = base.as_ref();
    let gof_path = append_suffix(path, ".gof");
    let file =
        File::open(&gof_path).with_context(|| format!("Failed to open GOF file {:?}", gof_path))?;
    let mut reader = BufReader::new(file);
    let mut entries = Vec::new();

    loop {
        match reader.read_u32::<LittleEndian>() {
            Ok(feature_id) => {
                let _pad = reader
                    .read_u32::<LittleEndian>()
                    .context("Failed reading GOF pad")?;
                let start = reader
                    .read_u64::<LittleEndian>()
                    .context("Failed reading GOF start offset")?;
                let end = reader
                    .read_u64::<LittleEndian>()
                    .context("Failed reading GOF end offset")?;
                entries.push(GofEntry {
                    feature_id,
                    start_offset: start,
                    end_offset: end,
                });
            }
            Err(e) if e.kind() == ErrorKind::UnexpectedEof => {
                break;
            }
            Err(e) => {
                return Err(e).context("Error reading GOF file");
            }
        }
    }

    Ok(entries)
}

pub fn load_gbi<P: AsRef<Path>>(base: P) -> Result<Vec<GbiEntry>> {
    let path = base.as_ref();
    let gbi_path = append_suffix(path, ".gbi");
    let mmap = safe_mmap_readonly(&gbi_path)?;
    let buf = &mmap[..];

    let mut cursor = 0;
    let seq_count = LittleEndian::read_u32(&buf[cursor..cursor + 4]) as usize;
    cursor += 4;

    let mut offsets = Vec::with_capacity(seq_count);
    for _ in 0..seq_count {
        let off = LittleEndian::read_u64(&buf[cursor..cursor + 8]) as usize;
        offsets.push(off);
        cursor += 8;
    }

    let mut entries = Vec::new();
    for &off in &offsets {
        let mut pos = off;
        let sid = LittleEndian::read_u32(&buf[pos..pos + 4]);
        pos += 4;

        let bin_count = LittleEndian::read_u32(&buf[pos..pos + 4]) as usize;
        pos += 4;

        for _ in 0..bin_count {
            let bin = LittleEndian::read_u32(&buf[pos..pos + 4]);
            pos += 4;
            let entry_count = LittleEndian::read_u32(&buf[pos..pos + 4]) as usize;
            pos += 4;

            for _ in 0..entry_count {
                let fid = LittleEndian::read_u32(&buf[pos..pos + 4]);
                pos += 4;
                let start = LittleEndian::read_u32(&buf[pos..pos + 4]);
                pos += 4;
                let end = LittleEndian::read_u32(&buf[pos..pos + 4]);
                pos += 4;

                entries.push(GbiEntry {
                    seqid_id: sid,
                    bin,
                    feature_id: fid,
                    start,
                    end,
                });
            }
        }
    }

    Ok(entries)
}

pub fn load_fts<P: AsRef<Path>>(gff_path: P) -> Result<Vec<String>> {
    let path = gff_path.as_ref();
    let fts_path = append_suffix(path, ".fts");
    let mmap = safe_mmap_readonly(&fts_path)?;
    let data = &mmap[..];

    let mut lines = Vec::new();
    let mut start = 0;

    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            let slice = &data[start..i];
            if !slice.is_empty() {
                let s = std::str::from_utf8(slice).context("FTS contains invalid UTF-8")?;
                lines.push(s.to_string());
            }
            start = i + 1;
        }
    }

    Ok(lines)
}

pub fn load_atn(path: &Path) -> Result<(String, Vec<String>)> {
    let atn_path = append_suffix(path, ".atn");
    let mmap = safe_mmap_readonly(&atn_path)?;
    let data = &mmap[..];
    let mut values = Vec::new();
    let mut attr_name = String::new();
    let mut start = 0;

    for (i, &b) in data.iter().enumerate() {
        if b == b'\n' {
            let slice = &data[start..i];
            if !slice.is_empty() {
                let line = std::str::from_utf8(slice)
                    .context("ATN contains invalid UTF-8")?
                    .trim();
                if let Some(stripped) = line.strip_prefix("#attribute=") {
                    attr_name = stripped.to_string();
                } else if !line.starts_with('#') {
                    values.push(line.to_string());
                }
            }
            start = i + 1;
        }
    }

    if attr_name.is_empty() {
        bail!("Missing #attribute=... header in .atn file");
    }

    Ok((attr_name, values))
}

pub fn load_prt<P: AsRef<Path>>(path: P) -> Result<Vec<PrtEntry>> {
    let path = path.as_ref();
    let prt_path = append_suffix(path, ".prt");
    let mmap = safe_mmap_readonly(&prt_path)?;
    if mmap.len() % 4 != 0 {
        bail!("Corrupted PRT: not aligned to u32");
    }

    let mut result = Vec::with_capacity(mmap.len() / 4);
    for (i, chunk) in mmap.chunks_exact(4).enumerate() {
        let parent = LittleEndian::read_u32(chunk);
        result.push(PrtEntry {
            child: i as u32,
            parent,
        });
    }
    Ok(result)
}

pub fn load_a2f<P: AsRef<Path>>(path: P) -> Result<Vec<A2fEntry>> {
    let path = path.as_ref();
    let a2f_path = append_suffix(path, ".a2f");
    let mmap = safe_mmap_readonly(&a2f_path)?;
    if mmap.len() % 4 != 0 {
        bail!("Corrupted A2F: not aligned to u32");
    }

    let mut result = Vec::with_capacity(mmap.len() / 4);
    for (i, chunk) in mmap.chunks_exact(4).enumerate() {
        let attr_id = LittleEndian::read_u32(chunk);
        result.push(A2fEntry {
            feature_id: i as u32,
            attr_id,
        });
    }
    Ok(result)
}
