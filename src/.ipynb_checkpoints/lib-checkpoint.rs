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
pub use index_loader::{
    core::{load_atn, load_sqs, safe_mmap_readonly},
    gof::{GofMap, load_gof},
    fts::{FtsMap, load_fts},
    prt::{PrtMap, load_prt},
    a2f::{A2fMap, load_a2f}
};


pub use utils::common::{
    CommonArgs, append_suffix, check_index_files_exist, write_gff_output, write_gff_output_filtered,
};
pub use utils::serial_interval_trees::{
    Interval, IntervalTree, save_multiple_trees, write_offsets_to_file,
};
