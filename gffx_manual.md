# GFFx Manual

GFFx is a high-performance Rust-based toolkit for fast extraction of features from large GFF3 annotation files. It supports both ID-based and region-based queries, with multi-threaded acceleration and a modular index system.

## Table of Contents

gffx <SUBCOMMAND> [OPTIONS]

Available subcommands:

- [index] Build index files
- [intersect] Extract features by region
- [extract] Extract features by ID
- [search] Search features by attribute
---

## Installation

cargo install gffx

Requires Rust >= 1.70.
---

## Subcommands

### index

Create index files from a GFF3 file.

```bash
gffx index -g genome.gff3
```

This generates all required `.gbi`, `.gof`, `.fts`, `.prt`, `.a2f`, `.atn`, `.sqs` files.

### extract

Extract feature models by feature ID(s).

```bash
gffx extract -g genome.gff3 --feature-id gene123
```

Options:
- `--feature-id <ID>`: A single feature ID.
- `--feature-file <FILE>`: A file containing multiple IDs (one per line).

### intersect

Extract feature models overlapping a region or regions.

```bash
gffx intersect -g genome.gff3 --bed regions.bed
```

Options:
- `--region chr1:1000-2000`: A single region.
- `--bed <FILE>`: A BED file containing regions.
- `--contained`, `--contains-region`, `--overlap`: Match mode.

### search

Search features by attribute values.

```bash
gffx search -g genome.gff3 --attr ID=gene123
```

Options:
- `--attr <KEY=VALUE>`: Search for a feature with a specific attribute.
- `--parent-of <ID>`: Find parent or ancestor models of a given feature.

---

## Index File Types

| File Extension | Purpose                            |
|----------------|------------------------------------|
| `.gbi`         | Bin-based region index             |
| `.gof`         | Byte offset index for GFF lines    |
| `.fts`         | Feature to sub-feature links       |
| `.prt`         | Parent-to-root model mapping       |
| `.a2f`         | Attribute to feature ID lookup     |
| `.atn`         | Attribute name table               |
| `.sqs`         | Sequence offset index              |

---

## License

GFFx is released under the MIT License.