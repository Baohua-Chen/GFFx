use anyhow::Result;
use clap::Parser;
use rayon::prelude::*;
use rand::seq::{IndexedRandom};
use rand::rng;
use std::{
    path::PathBuf,
};
use crate::{load_gof, write_gff_output};

/// Arguments
#[derive(Parser, Debug)]
#[command(
    about = "Sample feature groups per chromosome",
    long_about = "Sample feature groups per chromosome."
)]
pub struct SampleArgs {
    /// GFF file path (indexed via GOF)
    #[arg(short = 'i', long = "input", value_name = "FILE")]
    pub input: PathBuf,

    /// Ratio of downsampling
    #[arg(short = 'r', long = "ratio")]
    pub ratio: f32,
    
    /// Output file (required)
    #[arg(short = 'o', long = "output", value_name = "FILE")]
    pub output: Option<PathBuf>,

    /// Number of threads
    #[arg(short = 't', long = "threads", default_value_t = 12, value_name = "NUM")]
    pub threads: usize,

    /// Verbose logs
    #[arg(short = 'v', long = "verbose", default_value_t = false, value_name = "BOOL")]
    pub verbose: bool,
}

pub fn run(args: &SampleArgs) -> Result<()> {
    let verbose = args.verbose;
    let threads = if args.threads == 0 {
        std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1)
    } else {
        args.threads
    };
    let _ = rayon::ThreadPoolBuilder::new().num_threads(threads).build_global();
    let gff_path = &args.input;
    
    let gof = load_gof(&gff_path)?;

    let blocks: Vec<(u32, u64, u64)> = gof.seqid_index
        .par_iter()
        .flat_map(|(_seqid_num, indices)| {
            let mut rng = rng();
    
            // 1. collect all fids for this chromosome
            let fids: Vec<u32> = indices
                .iter()
                .map(|&i| gof.entries[i].feature_id)
                .collect();
    
            if fids.is_empty() {
                return Vec::new();
            }
    
            // 2. sample 10% fids
            let sample_size = (fids.len() as f32 * args.ratio).ceil() as usize;
            let sampled: Vec<u32> = fids.choose_multiple(&mut rng, sample_size).cloned().collect();
    
            // 3. use index_cached() to get offsets
            let idx = gof.index_cached();
            sampled
                .into_iter()
                .filter_map(|fid| idx.get(&fid).map(|&(s, e)| (fid, s, e)))
                .collect::<Vec<_>>()
        })
        .collect();

    // Step 3: write sampled GFF blocks
    write_gff_output(gff_path, &blocks, &args.output, verbose)?;
    Ok(())
}



        