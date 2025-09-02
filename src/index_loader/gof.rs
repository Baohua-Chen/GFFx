use anyhow::{bail, Context, Result};
use byteorder::{ByteOrder, LittleEndian};
use rustc_hash::FxHashMap;
use std::{path::Path, sync::OnceLock};
use crate::{append_suffix, safe_mmap_readonly};

#[derive(Debug)]
pub struct GofEntry {
    pub feature_id: u32,
    pub start_offset: u64,
    pub end_offset: u64,
}

#[derive(Debug)]
pub struct GofMap {
    /// Flat list of (feature_id, start, end) tuples as read from .gof
    pub entries: Vec<GofEntry>,
    /// Lazy, thread-safe cache: feature_id -> (start, end)
    index_cache: OnceLock<FxHashMap<u32, (u64, u64)>>,
}

impl GofMap {
    /// Build a transient index (allocates on every call).
    /// Prefer `index_cached()` for hot paths.
    pub fn index(&self) -> FxHashMap<u32, (u64, u64)> {
        self.entries
            .iter()
            .map(|e| (e.feature_id, (e.start_offset, e.end_offset)))
            .collect()
    }

    /// Get (or build once) the cached index.
    #[inline]
    pub fn index_cached(&self) -> &FxHashMap<u32, (u64, u64)> {
        self.index_cache.get_or_init(|| self.index())
    }

    /// O(1) lookup using the cached index.
    #[inline]
    pub fn get(&self, fid: u32) -> Option<&(u64, u64)> {
        self.index_cached().get(&fid)
    }

    /// Map a list of root IDs to their (start, end) offsets using the cached index.
    /// Missing roots are silently skipped (consistent with the original free function).
    #[inline]
    pub fn roots_to_offsets(
        &self,
        roots: &[u32],
        threads: usize,
    ) -> Vec<(u32, u64, u64)> {
        let idx = self.index_cached();

        // Simple heuristic: parallelize only for large inputs
        let should_parallel = threads > 1 && roots.len() > 2048;

        if should_parallel {
            use rayon::prelude::*;
            roots
                .par_iter()
                .filter_map(|&r| idx.get(&r).map(|&(s, e)| (r, s, e)))
                .collect()
        } else {
            let mut out = Vec::with_capacity(roots.len());
            for &r in roots {
                if let Some(&(s, e)) = idx.get(&r) {
                    out.push((r, s, e));
                }
            }
            out
        }
    }
}

/// Load a `.gof` file containing (u32 fid, u32 padding, u64 start, u64 end) records.
pub fn load_gof<P: AsRef<Path>>(gff_path: P) -> Result<GofMap> {
    let path = gff_path.as_ref();
    let gof_path = append_suffix(path, ".gof");
    let mmap = safe_mmap_readonly(&gof_path)
        .with_context(|| format!("Failed to mmap {}", gof_path.display()))?;
    let bytes = &mmap[..];
    const REC_SIZE: usize = 4 + 4 + 8 + 8;

    if bytes.len() % REC_SIZE != 0 {
        bail!(
            "Corrupted GOF ({}): length {} not multiple of {}",
            gof_path.display(),
            bytes.len(),
            REC_SIZE
        );
    }

    let mut entries = Vec::with_capacity(bytes.len() / REC_SIZE);
    for rec in bytes.chunks_exact(REC_SIZE) {
        let fid   = LittleEndian::read_u32(&rec[0..4]);
        // let _pad = LittleEndian::read_u32(&rec[4..8]);
        let start = LittleEndian::read_u64(&rec[8..16]);
        let end   = LittleEndian::read_u64(&rec[16..24]);
        entries.push(GofEntry { feature_id: fid, start_offset: start, end_offset: end });
    }

    Ok(GofMap {
        entries,
        index_cache: OnceLock::new(),
    })
}
