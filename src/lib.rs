// src/lib.rs
pub mod index_loader;
pub mod index_builder;
pub mod utils;
pub mod index;
//pub mod sort;
pub mod intersect;
pub mod extract;
pub mod search;

pub use index_builder::index_builder::build_index;
pub use index_loader::index_loader::{safe_mmap_readonly, GofEntry, PrtEntry, A2fEntry, load_gof, load_sqs, load_prt, load_fts, load_atn, load_a2f};
pub use utils::{
    utils::{CommonArgs, append_suffix, write_gff_output, check_index_files_exist},
    serial_interval_trees::{save_multiple_trees, write_offsets_to_file, IntervalTree, Interval},
};
