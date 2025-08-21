// src/lib.rs
pub mod index;
pub mod index_builder;
pub mod index_loader;
pub mod utils;
//pub mod sort;
pub mod extract;
pub mod intersect;
pub mod search;

pub use index_builder::core::build_index;
pub use index_loader::core::{
    A2fEntry, GofEntry, PrtEntry, load_a2f, load_atn, load_fts, load_gof, load_prt, load_sqs,
    safe_mmap_readonly,
};
pub use utils::common::{CommonArgs, append_suffix, check_index_files_exist, write_gff_output};
pub use utils::serial_interval_trees::{
    Interval, IntervalTree, save_multiple_trees, write_offsets_to_file,
};
