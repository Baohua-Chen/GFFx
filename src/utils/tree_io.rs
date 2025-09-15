use crate::utils::tree::IntervalTree;
use anyhow::{bail, Context, Result};
use bincode2::{deserialize, deserialize_from, serialize, serialize_into};

use memmap2::MmapOptions;
use serde::de::DeserializeOwned;
use std::{
    fs::{self, File},
    io::{BufReader, Read, Seek, SeekFrom, BufWriter, Write},
    path::Path,
};

impl<T> IntervalTree<T>
where
    T: Ord + Copy + serde::Serialize + for<'de> serde::Deserialize<'de>,
{
    /// Serialize the whole tree to a file via bincode2.
    pub fn save_to_file(&self, path: &Path) -> std::io::Result<()> {
        let encoded = serialize(self).map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        fs::write(path, encoded)?;
        Ok(())
    }

    /// Deserialize a tree from a file (whole-file read).
    pub fn load_from_file(path: &Path) -> std::io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buf = Vec::new();
        file.read_to_end(&mut buf)?;
        let tree: Self = deserialize(&buf)
            .map_err(|e| std::io::Error::new(std::io::ErrorKind::Other, e))?;
        Ok(tree)
    }
}

/// Save multiple trees back-to-back into a single `.rit` file.
/// Offsets of each tree (in bytes) should be recorded separately.
pub fn save_multiple_trees<T>(trees: &[IntervalTree<T>], rit_path: &Path) -> Result<Vec<u64>>
where
    T: Ord + Copy + serde::Serialize + for<'de> serde::Deserialize<'de>,
{
    let mut offsets = Vec::with_capacity(trees.len());
    let file = File::create(rit_path)?;
    let mut writer = BufWriter::new(file);

    for tree in trees {
        // Record current offset
        let pos = writer.seek(SeekFrom::Current(0))?;
        offsets.push(pos);
        // Write tree
        serialize_into(&mut writer, tree)?;
    }
    writer.flush()?;
    Ok(offsets)
}

/// Write offsets (Vec<u64>) as JSON into `.rix` file.
pub fn write_offsets_to_file(offsets: &[u64], rix_path: &Path) -> Result<()> {
    let file = File::create(rix_path)?;
    let mut writer = BufWriter::new(file);
    serde_json::to_writer(&mut writer, offsets)?;
    writer.flush()?;
    Ok(())
}

/// Load multiple trees with streaming (using BufReader).
pub fn load_trees_streaming<T>(rit_path: &Path, rix_path: &Path) -> Result<Vec<IntervalTree<T>>>
where
    T: Ord + Copy + DeserializeOwned,
{
    let mut reader = BufReader::new(File::open(rit_path)?);
    let offsets: Vec<u64> = {
        let f = File::open(rix_path)?;
        serde_json::from_reader(f)?
    };

    let mut trees = Vec::with_capacity(offsets.len());
    for &off in &offsets {
        reader.seek(SeekFrom::Start(off))?;
        let tree: IntervalTree<T> = deserialize_from(&mut reader)
            .map_err(|e| anyhow::anyhow!("Deserializing failed: {}", e))?;
        trees.push(tree);
    }
    Ok(trees)
}

/// Load multiple trees via memory-mapping.
pub fn load_trees_mmap<T>(rit_path: &Path, rix_path: &Path) -> Result<Vec<IntervalTree<T>>>
where
    T: Ord + Copy + for<'de> serde::Deserialize<'de>,
{
    let file = File::open(rit_path).with_context(|| format!("open {}", rit_path.display()))?;
    let mmap = unsafe { MmapOptions::new().map(&file) }
        .with_context(|| format!("mmap {}", rit_path.display()))?;
    let buf: &[u8] = &mmap;

    let offsets: Vec<u64> = {
        let f = File::open(rix_path).with_context(|| format!("open {}", rix_path.display()))?;
        serde_json::from_reader::<_, Vec<u64>>(f)
            .with_context(|| format!("parse json {}", rix_path.display()))?
    };

    if offsets.is_empty() {
        return Ok(Vec::new());
    }

    for w in offsets.windows(2) {
        if w[0] > w[1] {
            bail!("offsets not sorted ascending: {:?} > {:?}", w[0], w[1]);
        }
    }
    let last = *offsets.last().unwrap() as usize;
    if last > buf.len() {
        bail!("last offset {} out of file size {}", last, buf.len());
    }

    let mut out = Vec::with_capacity(offsets.len());
    for (i, start_u64) in offsets.iter().copied().enumerate() {
        let start = start_u64 as usize;
        let end = if i + 1 < offsets.len() {
            offsets[i + 1] as usize
        } else {
            buf.len()
        };
        if end < start || end > buf.len() {
            bail!("bad slice range: {}..{} (file len {})", start, end, buf.len());
        }
        let slice = &buf[start..end];
        let tree: IntervalTree<T> = deserialize(slice)
            .with_context(|| format!("bincode2 deserialize tree #{} ({}..{})", i, start, end))?;
        out.push(tree);
    }

    Ok(out)
}
