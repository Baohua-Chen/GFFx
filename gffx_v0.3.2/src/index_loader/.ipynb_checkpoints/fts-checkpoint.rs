    use anyhow::{Context, Result};
    use rustc_hash::{FxHashMap, FxHashSet};
    use std::path::Path;
    use std::sync::OnceLock; // for lazy cached index
    use rayon::prelude::*;    // only needed if you want the parallel path inside methods
    use crate::{append_suffix, safe_mmap_readonly};
    
    #[derive(Debug)]
    pub struct FtsMap {
        pub ids: Vec<String>,
        /// Lazy-built cache: String -> numeric ID (u32)
        ///
        /// We cache an owned-key map to avoid lifetime issues with &str keys.
        index_cached: OnceLock<FxHashMap<String, u32>>,
    }
    
    impl FtsMap {
        /// Build an owned-key index (String -> u32). Used by the cache.
        fn build_owned_index(&self) -> FxHashMap<String, u32> {
            let mut m = FxHashMap::with_capacity_and_hasher(self.ids.len(), Default::default());
            for (i, s) in self.ids.iter().enumerate() {
                m.insert(s.clone(), i as u32);
            }
            m
        }
    
        /// Get the cached index (String -> u32), building it once if necessary.
        pub fn index(&self) -> &FxHashMap<String, u32> {
            self.index_cached.get_or_init(|| self.build_owned_index())
        }
    
        /// Return numeric ID for a given string ID (via cached index).
        pub fn get_fid(&self, id: &str) -> Option<u32> {
            self.index().get(id).copied()
        }
    
        /// Return string ID for a given numeric ID.
        pub fn get_id(&self, fid: u32) -> Option<&str> {
            self.ids.get(fid as usize).map(|s| s.as_str())
        }
    
        /// Map a set of string feature IDs to numeric IDs, returning:
        /// - a set of numeric IDs that were found
        /// - a vector of missing string IDs (not present in the index)
        ///
        /// It selects a parallel or sequential path based on `threads` and input size.
        pub fn map_fnames_to_fids(
            &self,
            feature_names: &FxHashSet<String>,
            threads: usize,
        ) -> (FxHashSet<u32>, Vec<String>) {
            // Small enum to separate successful mappings from missing IDs
            enum Either<L, R> { Left(L), Right(R) }
    
            let idx = self.index();
    
            // Closure: convert each string ID into Either::Left(fid) or Either::Right(fname)
            let mapper = |fname: &String| {
                if let Some(&fid) = idx.get(fname) {
                    Either::Left(fid)
                } else {
                    Either::Right(fname.clone())
                }
            };
    
            if threads > 1 && feature_names.len() > 2 {
                // Parallel execution (Rayon)
                feature_names
                    .par_iter()
                    .map(mapper)
                    .fold(
                        || (FxHashSet::default(), Vec::new()),
                        |mut acc, e| {
                            match e {
                                Either::Left(n)  => { acc.0.insert(n); } // collect numeric ID
                                Either::Right(s) => { acc.1.push(s); }   // collect missing string
                            }
                            acc
                        },
                    )
                    .reduce(
                        || (FxHashSet::default(), Vec::new()),
                        |mut a, b| {
                            a.0.extend(b.0);
                            a.1.extend(b.1);
                            a
                        },
                    )
            } else {
                // Sequential execution
                let mut set = FxHashSet::default();
                set.reserve(feature_names.len());
                let mut miss = Vec::new();
                for fname in feature_names {
                    if let Some(&fid) = idx.get(fname) {
                        set.insert(fid);
                    } else {
                        miss.push(fname.clone());
                    }
                }
                (set, miss)
            }
        }
    }
    
    /// Load `.fts` file into FtsMap (unchanged except initializing the cache).
    pub fn load_fts<P: AsRef<Path>>(gff_path: P) -> Result<FtsMap> {
        let path = gff_path.as_ref();
        let fts_path = append_suffix(path, ".fts");
    
        let mmap = safe_mmap_readonly(&fts_path)
            .with_context(|| format!("Failed to mmap {}", fts_path.display()))?;
        let data = &mmap[..];
    
        let mut lines = Vec::new();
        let mut start = 0;
    
        for (i, &b) in data.iter().enumerate() {
            if b == b'\n' {
                let mut slice = &data[start..i];
                if !slice.is_empty() {
                    if slice.ends_with(b"\r") {
                        slice = &slice[..slice.len() - 1];
                    }
                    let s = std::str::from_utf8(slice)
                        .with_context(|| format!("FTS contains invalid UTF-8 at byte {}", start))?;
                    lines.push(s.to_string());
                }
                start = i + 1;
            }
        }
        if start < data.len() {
            let mut slice = &data[start..];
            if !slice.is_empty() {
                if slice.ends_with(b"\r") {
                    slice = &slice[..slice.len() - 1];
                }
                let s = std::str::from_utf8(slice)
                    .with_context(|| format!("FTS contains invalid UTF-8 at byte {}", start))?;
                lines.push(s.to_string());
            }
        }
    
        Ok(FtsMap { ids: lines, index_cached: OnceLock::new() })
    }