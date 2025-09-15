use crate::{IntervalTree, load_sqs, append_suffix};
use anyhow::{bail, Context, Result};
use bincode2::deserialize;
use memmap2::MmapOptions;
use rustc_hash::FxHashMap;
use std::{fs::File, path::Path};

/// Application-facing structure:
/// - per-sequence interval trees
/// - string -> numeric ID mapping
#[derive(Debug)]
pub struct TreeIndexData {
    pub chr_entries: FxHashMap<u32, IntervalTree<u32>>,
    pub seqid_to_num: FxHashMap<String, u32>,
    pub num_to_seqid: Vec<String>,
}

impl TreeIndexData {
    /// Construct TreeIndexData from explicit sequence IDs and `.rit/.rix` files.
    /// Convenience: derive `{prefix}.rit/.rix` from `index_prefix`.
    pub fn load_tree_index<P: AsRef<Path>>(gff_path: P) -> Result<Self> {
        let path = gff_path.as_ref();
        let (num_to_seqid, seqid_to_num) = load_sqs(path)?;
        let rit_path = append_suffix(path, ".rit");
        let rix_path = append_suffix(path, ".rix");
        
        let chr_entries = Self::load_region_index(&rit_path, &rix_path)?;
        
        Ok(Self {
            chr_entries,
            seqid_to_num,
            num_to_seqid
        })
    }

    fn load_region_index(
        rit_path: &Path,
        rix_path: &Path,
    ) -> Result<FxHashMap<u32, IntervalTree<u32>>> {
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
            return Ok(FxHashMap::default());
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

        let mut map = FxHashMap::with_capacity_and_hasher(offsets.len(), Default::default());
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
            let tree: IntervalTree<u32> = deserialize(slice).with_context(|| {
                format!("bincode2 deserialize tree #{} ({}..{})", i, start, end)
            })?;
            map.insert(i as u32, tree);
        }
        Ok(map)
    }
}
