use anyhow::{bail, Result};
use byteorder::{ByteOrder, LittleEndian};
use rustc_hash::FxHashMap;
use std::{path::Path, sync::OnceLock};
use rayon::prelude::*; // Parallel iteration (no feature gate)

use crate::{append_suffix, safe_mmap_readonly};

#[derive(Debug, Clone, Copy)]
pub struct PrtEntry {
    /// Child node id
    pub child: u32,
    /// Parent node id
    pub parent: u32,
}

#[derive(Debug)]
pub struct PrtMap {
    /// Flat list of (child -> parent) entries, where index i is the child id
    pub entries: Vec<PrtEntry>,
    /// Lazy, thread-safe cache for child->parent index
    index_cache: OnceLock<FxHashMap<u32, u32>>,
}

impl PrtMap {
    /// Construct a PrtMap from a list of entries.
    pub fn new(entries: Vec<PrtEntry>) -> Self {
        Self {
            entries,
            index_cache: OnceLock::new(),
        }
    }

    /// Build a child -> parent hashmap for O(1) lookups (allocates every call).
    pub fn index(&self) -> FxHashMap<u32, u32> {
        self.entries.iter().map(|e| (e.child, e.parent)).collect()
    }

    /// Get the cached child -> parent index, building it once on first use.
    pub fn index_cached(&self) -> &FxHashMap<u32, u32> {
        self.index_cache.get_or_init(|| self.index())
    }

    /// Small helper: get the parent of a child via linear scan (O(n)).
    pub fn get_parent(&self, child: u32) -> Option<u32> {
        self.entries
            .iter()
            .find(|e| e.child == child)
            .map(|e| e.parent)
    }

    /// Resolve the root of a node by following parent pointers via *array access* (fast path).
    #[inline]
    fn resolve_root(&self, start: u32) -> (u32, bool) {
        let n = self.entries.len() as u32;
        let mut cur = start;

        loop {
            if cur >= n {
                return (u32::MAX, true);
            }
            let p = self.entries[cur as usize].parent;

            if p == cur {
                return (cur, false);
            }
            if p >= n {
                return (u32::MAX, true);
            }
            cur = p;
        }
    }

    #[inline]
    /// Map a Vec<FID> to a Vec<ROOT> using the fast resolver.
    /// - If a FID equals `u32::MAX`, keep it as-is (sentinel).
    /// - Output order matches the input order.
    pub fn map_fids_to_roots(&self, fids: &Vec<u32>, threads: usize) -> Vec<u32> {
        let should_parallel = threads > 1 && fids.len() > 256; // tune threshold as needed

        if should_parallel {
            fids.par_iter()
                .map(|&fid| {
                    if fid == u32::MAX {
                        u32::MAX
                    } else {
                        self.resolve_root(fid).0
                    }
                })
                .collect()
        } else {
            let mut out = Vec::with_capacity(fids.len());
            for &fid in fids {
                if fid == u32::MAX {
                    out.push(u32::MAX);
                } else {
                    out.push(self.resolve_root(fid).0);
                }
            }
            out
        }
    }
}

/// Load a `.prt` file that encodes parent pointers as a u32 array.
/// Each 4-byte little-endian word is the parent id of the child at the same index.
/// For child i, parent = data[i].
pub fn load_prt<P: AsRef<Path>>(gff_path: P) -> Result<PrtMap> {
    let path = gff_path.as_ref();
    let prt_path = append_suffix(path, ".prt");
    let mmap = safe_mmap_readonly(&prt_path)?;
    if mmap.len() % 4 != 0 {
        bail!("Corrupted PRT: not aligned to u32");
    }

    let mut entries = Vec::with_capacity(mmap.len() / 4);
    for (i, chunk) in mmap.chunks_exact(4).enumerate() {
        let parent = LittleEndian::read_u32(chunk);
        entries.push(PrtEntry {
            child: i as u32,
            parent,
        });
    }
    Ok(PrtMap::new(entries))
}
