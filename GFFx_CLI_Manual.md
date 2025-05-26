
# GFFx Command Line Manual

**GFFx** is a high-performance, Rust-based toolkit for extracting and querying annotations from GFF3 files. It supports fast indexing and feature retrieval with several subcommands.

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

## `index`: Build Index for GFF

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

## `intersect`: Extract Features by Region

Extracts models intersecting with regions from a GFF file, either from a single region or a BED file.

```bash
gffx intersect [OPTIONS] --input <INPUT> <--region <REGION>|--bed <BED>>
```

**Options:**
Required
| Option                      | Description                                                  |
| --------------------------- | ------------------------------------------------------------ |
| `-i`, `--input` `<INPUT>`   | Input GFF file path                                          |
| *(one of)*                  |                                                              |
| `-r`, `--region` `<REGION>` | Single region in `chr:start-end` format                      |
| `-b`, `--bed` `<BED>`       | BED file containing multiple regions                         |

Optional
| Option                      | Description                                                  |
| -------------------------   | ------------------------------------------------------------ |
| `-o`, `--out` `<OUT>`       | Output file path (default: stdout)                           |
| `-T`, `--types` `<TYPES>`   | Comma-separated feature types to retain (e.g., `gene,exon`)  |
| `-t`, `--threads` `<N>`     | Number of threads to use (default: 4)                        |
| `-c`, `--contained`         | Only keep features fully contained within the region         |
| `-C`, `--contains-region`   | Only keep features that fully contain the region             |
| `-O`, `--overlap`           | Keep features that partially or fully overlap (default mode) |
| `-v`, `--invert`            | Invert selection (exclude matched features)                  |
| `-V`, `--verbose`           | Enable verbose output                                        |
| `-h`, `--help`              | Show help message                                            |

---

## `extract`: Extract Features by ID

Extracts annotation models by feature ID(s), including their parent models.

```bash
gffx extract --feature-file ids.txt -i input.gff
```

**Options:**

| Option                    | Description                              |
|---------------------------|------------------------------------------|
| `--feature-id`            | Single feature ID                        |
| `--feature-file`          | File containing feature IDs              |
| `-i`, `--input`           | Input GFF file                           |
| `-o`, `--output`          | Output GFF file                          |

---

## 🔍 `search`: Search Features by Attribute

Searches for features using a specified attribute value and retrieves the full annotation models.

```bash
gffx search -a geneX -i input.gff
```

**Options:**

| Option                    | Description                              |
|---------------------------|------------------------------------------|
| `--attr`, `-a`            | Attribute value to search for            |
| `--attr-file`             | File containing attribute values         |
| `--threads`               | Number of threads to use                 |
| `-i`, `--input`           | Input GFF file                           |
| `-o`, `--output`          | Output GFF file                          |

---

## 💡 Example Use Cases

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

## 🛠️ Notes

- Make sure you run `gffx index` before using `intersect`, `extract`, or `search`.
- All subcommands are optimized for multi-threaded execution where applicable.
