
# GFFx Command Line Manual

**GFFx** is a high-performance, Rust-based toolkit for extracting and querying annotations from GFF3 files. It supports fast indexing and feature retrieval with several subcommands.
It can be used both as a **command-line tool** and as a **Rust library**.

## Table of Contents

*GFFx version 0.2.0a*

- [Installation](#installation)
- [Basic Usage](#basic-usage)
  - [index](#index) - Build index files
  - [extract](#extract) - Extract features by ID
  - [intersect](#intersect) - Extract features by regions
  - [search](#search) - Search features by attributes
- [Example Use Cases](#example-use-cases)
- [Using GFFx as a Rust Library](#using-gffx-as-a-rust-library)
- [Index File Types](#index-file-types)
- [License](#license)
- [Notes](#notes)

---
## Installation

### Option 1: Install via [crates.io](https://crates.io/crates/gffx)

```bash
cargo install gffx                  # install to default location (~/.cargo/bin)
cargo install gffx --root /your/path  # optional: install to custom location
```

### Option 2: Install from source

```bash
git clone https://github.com/Baohua-Chen/GFFx.git
cd GFFx
cargo build --release
# Optionally copy the binary
cp target/release/gffx /your/path
```

> Requires **Rust 1.70 or later**. You can install or update Rust using [rustup](https://rustup.rs).
---


## Basic Usage

```bash
gffx <SUBCOMMAND> [OPTIONS]
```

Available subcommands:

- [index] Build index files
- [intersect] Extract features by region
- [extract] Extract features by ID
- [search] Search features by attribute

---

### `index`

Builds index files from a GFF file to accelerate downstream operations.

```bash
gffx index [OPTIONS] --input <INPUT>
```

**Options:**

| Option                 | Description                                     |
|------------------------|-------------------------------------------------|
| `-i`, `--input`        | Input GFF file                                  |
| `-a`, `--attribute`    | Attribute key to extract (default: `gene_name`) |
| `-v`, `--verbose`      | Enable verbose output                           |
| `-h`, `--help   `      | Print help                                      |

---

### `intersect`

Extracts models intersecting with regions from a GFF file, either from a single region or a BED file.

```bash
gffx intersect [OPTIONS] --input <INPUT> <--region <REGION>|--bed <BED>>
```

**Options:**
Required
| Option                      | Description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| `-i`, `--input` `<INPUT>`   | Input GFF file path                                          |
| `-r`, `--region` `<REGION>` | Single region in `chr:start-end` format                      |
| `-b`, `--bed` `<BED>`       | BED file containing multiple regions                         |

> **Note**: Exactly one of `--region` or `--bed` must be specified.


Optional
| Option                      | Description                                                                    |
| -------------------------   | ------------------------------------------------------------------------------ |
| `-o`, `--output` `<OUT>`       | Output file path (default: stdout)                                          |
| `-v`, `--invert`            | Invert selection (exclude matched features)                                    |
| `-T`, `--types` `<TYPES>`   | Filter output to include only features of specified types (e.g., `gene,exon`)  |
| `-V`, `--verbose`           | Enable verbose output                                                          |
| `-h`, `--help`              | Show help message                                                              |
| *(one of)*                  |                                                                                |
| `-c`, `--contained`         | Only keep features fully contained within the region                           |
| `-C`, `--contains-region`   | Only keep features that fully contain the region                               |
| `-O`, `--overlap`           | Keep features that partially or fully overlap (default mode)                   |

---

### `extract`

Extracts annotation models by feature ID(s), including their parent models.

```bash
gffx extract [OPTIONS] --input <INPUT> <--feature-file <FEATURE_FILE>|--feature-id <FEATURE_ID>>```
```

**Options:**

Required
| Option                                   | Description                                                  |
| ---------------------------------------- | ------------------------------------------------------------ |
| `-i`, `--input` `<INPUT>`                | Input GFF file path                                          |
| *(one of)*                               |                                                              |
| `-f`, `--feature-id` `<FEATURE_ID>`      | Extrach by a single feature id                               |
| `-F`, `--feature-file` `<FEATURE_FILE>`  | Extrach by a BED file containing multiple regions            |

Optional
| Option                      | Description                                                                    |
| -------------------------   | ------------------------------------------------------------------------------ |
| `-o`, `--output` `<OUT>`    | Output file path (default: stdout)                                             |
| `-T`, `--types` `<TYPES>`   | Filter output to include only features of specified types (e.g., `gene,exon`)  |
| `-d`, `--descendants-only`  | Only extract feature(s) and their/its descendants                              |
| `-V`, `--verbose`           | Enable verbose output                                                          |
| `-h`, `--help`              | Show help message                                                              |
---

### `search`

Searches for features using a specified attribute value and retrieves the full annotation models.

```bash
gffx search -a geneX -i input.gff
```

**Options:**

Required
| Option                                   | Description                                                  |
| ---------------------------------------- | ------------------------------------------------------------ |
| `-i`, `--input` `<INPUT>`                | Input GFF file path                                          |
| *(one of)*                               |                                                              |
| `-a`, `--attr` `ATTRIBUTE_VALUE>`        | Search a single attribute value/pattern                      |
| `-A`, `--attr-list` `<ATTRIBUTE_LIST>`   | Search attribute values/patterns defined in a text file      |

Optional
| Option                      | Description                                                                    |
| -------------------------   | ------------------------------------------------------------------------------ |
| `-o`, `--output` `<OUT>`    | Output file path (default: stdout)                                             |
| `-r`, `--regex` `<REGEX>`   | Enable regex matching for attribute values                                     |
| `-T`, `--types` `<TYPES>`   | Filter output to include only features of specified types (e.g., `gene,exon`)  |
| `-V`, `--verbose`           | Enable verbose output                                                          |
| `-h`, `--help`              | Show help message                                                              |
---

## Example Use Cases

```bash
# Build index
gffx index -i genes.gff -a gene_name

# Extract all features overlapping with a region
gffx intersect --region chr1:10000-20000 -i genes.gff -o out.gff

# Extract models from a list of gene IDs
gffx extract --feature-file genes.txt -i genes.gff -o subset.gff

# Search by gene name and extract the full model
gffx search -a TP53 -i genes.gff -o tp53_model.gff
```

---

## Using GFFx as a Rust Library

You can use GFFx as a Rust library in your own project.

### Add to Cargo.toml

```toml
[dependencies]
gffx = "0.1.1"
```

### Example: Manually extract features from region using index files
The following example runs inside a main() -> Result<()> context:

```rust
use anyhow::Result;
use gffx::{
    CommonArgs, IntersectArgs, IndexData, parse_region, query_features,
    extract_gff_blocks, load_gof
};
use std::path::PathBuf;
use rustc_hash::FxHashMap;

fn main() -> Result<()> {
    let input_path = PathBuf::from("example.gffx");
    let region_str = "chr1:1000-2000";

    let common = CommonArgs {
        input: input_path.clone(),
        output: None, // Write to stdout
        verbose: true,
    };

    // Load sequence ID map
    let (_, seqid_map) = gffx::load_sqs(&input_path)?;

    // Parse the region string
    let region = parse_region(&region_str, &seqid_map, &common)?;

    // Load the index data (interval trees)
    let index_data = IndexData::load(&input_path, &common)?;

    // Query features overlapping the region
    let feats = query_features(&index_data, vec![region], false, false, false)?;

    // Collect unique feature IDs
    let mut ids: Vec<u32> = feats.iter().map(|&(id, _, _)| id).collect();
    ids.sort_unstable(); ids.dedup();

    // Load GFF Offset Format (.gof) file
    let gof = load_gof(&input_path)?;
    let gof_map: FxHashMap<_, _> = gof.into_iter()
        .map(|e| (e.feature_id, (e.start_offset, e.end_offset)))
        .collect();

    // Extract GFF blocks from the indexed GFFx file
    extract_gff_blocks(&input_path, &gof_map, &ids, &None, true)?;

    Ok(())
}
```

---

### Available Public APIs

#### Index building & checking (`index_builder`)
- `build_index`

#### Index loading (`index_loader`)
- `load_gof`, `load_prt`, `load_fts`, `load_atn`, `load_a2f`, `load_sqs`
- `safe_mmap_readonly`
- `GofEntry`, `PrtEntry`, `A2fEntry`

#### Interval querying data structures (`utils::serial_interval_trees`)
- `IntervalTree`, `Interval`
- `save_multiple_trees`, `write_offsets_to_file`

#### Other utilities (`utils`)
- `write_gff_output`
- `check_index_files_exist`
- `append_suffix`

#### Command-line compatibility
- `CommonArgs`

---


## Index File Types

| File Extension | Purpose                                             |
|----------------|-----------------------------------------------------|
| `.gof`         | Byte offset index for GFF feature blocks            |
| `.fts`         | Feature ID table                                    |
| `.prt`         | Child to parent mapping                             |
| `.a2f`         | Attribute to feature ID mapping                     |
| `.atn`         | Attribute value table                               |
| `.sqs`         | Sequence ID table                                   |
| `.rit`         | Interval tree index                                 |
| `.rix`         | Byte offest index for interval trees in.rit file    |

---

## Notes

- Make sure you run `gffx index` before using `intersect`, `extract`, or `search`.

---

## License

GFFx is released under the MIT or Apache-2.0 License.

---
