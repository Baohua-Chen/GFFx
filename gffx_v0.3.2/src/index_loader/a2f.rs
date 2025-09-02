use anyhow::{bail, Context, Result};
use byteorder::{ByteOrder, LittleEndian};
use rustc_hash::{FxHashMap, FxHashSet};
use std::path::Path;

use crate::{append_suffix, safe_mmap_readonly};

/// A2fMap stores two indexes:
/// - `aid_to_fids`: mapping from Attribute ID (AID) to all Feature IDs (FIDs) that reference it.
/// - `fid_to_aid`: reverse mapping from Feature ID (FID) to its Attribute ID (if any).
///
/// The `.a2f` file on disk is encoded as `fid -> aid` (one u32 per fid).
/// During loading, we construct the reverse map `aid -> fids` for efficient queries.
#[derive(Debug)]
pub struct A2fMap {
    aid_to_fids: FxHashMap<u32, Vec<u32>>,
    fid_to_aid: Vec<Option<u32>>,
}

impl A2fMap {
    pub fn new(aid_to_fids: FxHashMap<u32, Vec<u32>>, fid_to_aid: Vec<Option<u32>>) -> Self {
        Self { aid_to_fids, fid_to_aid }
    }

    /// Get all FIDs associated with a given AID.
    #[inline]
    pub fn fids_for_aid(&self, aid: u32) -> Option<&[u32]> {
        self.aid_to_fids.get(&aid).map(|v| v.as_slice())
    }

    /// Get the AID associated with a given FID (if any).
    #[inline]
    pub fn aid_for_fid(&self, fid: u32) -> Option<u32> {
        self.fid_to_aid.get(fid as usize).and_then(|x| *x)
    }

    /// Map a set of AIDs into a set of FIDs (deduplicated, no order guaranteed).
    #[inline]
    pub fn map_aids_to_fids_set(&self, aids: &FxHashSet<u32>) -> FxHashSet<u32> {
        let mut out = FxHashSet::default();
        for &aid in aids {
            if let Some(fids) = self.aid_to_fids.get(&aid) {
                out.extend(fids.iter().copied());
            } else {
                eprintln!("[WARN] AID {} not found (no FIDs).", aid);
            }
        }
        out
    }

    /// Map a list of AIDs into a combined Vec of FIDs.
    /// The caller should sort/deduplicate if stable order is needed.
    #[inline]
    pub fn map_aids_to_fids_vec(&self, aids: &[u32]) -> Vec<u32> {
        let mut out = Vec::new();
        for &aid in aids {
            if let Some(fids) = self.aid_to_fids.get(&aid) {
                out.extend_from_slice(fids);
            } else {
                eprintln!("[WARN] AID {} not found (no FIDs).", aid);
            }
        }
        out
    }

    #[inline]
    pub fn len_fids(&self) -> usize { self.fid_to_aid.len() }

    #[inline]
    pub fn is_empty(&self) -> bool { self.fid_to_aid.is_empty() }
}

/// Load `.a2f` file and build A2fMap:
/// - On disk: each 4-byte little-endian u32 represents the AID for a given FID.
/// - Value `u32::MAX` means "no attribute" (None).
/// - In memory: build both `fid -> aid` (vector) and `aid -> fids` (hashmap).
pub fn load_a2f<P: AsRef<Path>>(gff_path: P) -> Result<A2fMap> {
    let path = gff_path.as_ref();
    let a2f_path = append_suffix(path, ".a2f");

    let mmap = safe_mmap_readonly(&a2f_path)
        .with_context(|| format!("Failed to mmap {}", a2f_path.display()))?;

    if mmap.len() % 4 != 0 {
        bail!(
            "Corrupted A2F ({}): length {} not aligned to u32",
            a2f_path.display(),
            mmap.len()
        );
    }

    let n = mmap.len() / 4;
    let mut fid_to_aid: Vec<Option<u32>> = Vec::with_capacity(n);
    let mut aid_to_fids: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

    for (fid, chunk) in mmap.chunks_exact(4).enumerate() {
        let raw = LittleEndian::read_u32(chunk);
        let aid = if raw == u32::MAX { None } else { Some(raw) };
        fid_to_aid.push(aid);

        if let Some(aid_val) = aid {
            aid_to_fids.entry(aid_val).or_default().push(fid as u32);
        }
    }

    // Optional: deduplicate and sort FID lists
    for fids in aid_to_fids.values_mut() {
        fids.sort_unstable();
        fids.dedup();
    }

    Ok(A2fMap::new(aid_to_fids, fid_to_aid))
}
