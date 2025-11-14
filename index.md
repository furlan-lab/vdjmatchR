# vdjmatchR

![vdjmatchR logo](reference/figures/logo2.png)

**Fast TCR sequence matching with VDJdb in R**

[Documentation](https://furlan-lab.github.io/vdjmatchR/) •
[Installation](#installation) • [Quickstart](#quickstart) •
[Vignettes](#vignettes)

------------------------------------------------------------------------

## Overview

vdjmatchR provides a high-performance R interface to TCR sequence
matching against the VDJdb database. Built with a Rust core via
[`extendr`](https://extendr.github.io/), it combines the speed of Rust
with the convenience of R.

### Key Features

✨ **Fast Matching** - Rust-powered parallel processing for high
throughput 🔍 **Fuzzy Matching** - Configurable edit distance tolerance
📊 **Database Tools** - Easy filtering, exploration, and conversion to
data.table 🧬 **Seurat Integration** - Seamless integration with 10x
Genomics VDJ data 📈 **Progress Tracking** - Real-time progress bars for
large datasets 🎯 **VDJdb Compatible** - Works with both slim and full
VDJdb formats

## Installation

### Prerequisites

- R (\>= 4.0)
- Rust (\>= 1.70) - [Install
  Rust](https://www.rust-lang.org/tools/install)
- R development tools (Rtools on Windows, Xcode Command Line Tools on
  macOS)

### From Github in R

``` r
devtools::install_github("furlan-lab/vdjmatchR")
```

### Clone and install (from source) in shell

``` sh
git clone https://github.com/furlan-lab/vdjmatchR.git
cd vdjmatchR
R CMD INSTALL .
```

## Quickstart

``` r
library(vdjmatchR)

# Load the packaged VDJdb database
db_path <- vdjdb_packaged_path(use_fat_db = TRUE)
db <- vdjdb_open_file(db_path)

# Filter for human TRB entries with high confidence
db_filtered <- filter_db(db,
                         species = "HomoSapiens",
                         gene = "TRB",
                         min_vdjdb_score = 2)

# Match a single TCR (exact match)
result <- match_tcr_df(
  db_filtered,
  cdr3 = "CASSLGQAYEQYF",
  v_segment = "TRBV12-3",
  j_segment = "TRBJ2-7",
  scope = "0,0,0,0",  # exact match
  top_n = 5
)

# Batch matching with fuzzy tolerance (parallel processing enabled)
cdr3s <- c("CASSLGQAYEQYF", "CASSIRSSYEQYF", "CASSPGQGTTDTQYF")
vs    <- c("TRBV12-3", "TRBV7-9", "TRBV5-4")
js    <- c("TRBJ2-7", "TRBJ2-7", "TRBJ2-3")

results <- match_tcr_many_df(
  db_filtered,
  cdr3s, vs, js,
  scope = "3,1,2,3",  # fuzzy: max 3 substitutions, 1 insertion, 2 deletions
  top_n = 3,
  progress = TRUE  # shows progress bar
)

head(results)
```

## Core Functions

### Database Operations

- **`vdjdb_open_file(path)`** - Load VDJdb database from file
- **`vdjdb_packaged_path(use_fat_db)`** - Get path to bundled database
- **`filter_db(db, species, gene, min_vdjdb_score)`** - Filter database
  entries
- **`filter_db_by_epitope_size(db, min_size)`** - Filter by epitope
  representation
- **`db_to_table(db)`** - Convert database to data.table for exploration
- **`db_summary(db)`** - Print database summary statistics

### TCR Matching

- **`match_tcr_df(db, cdr3, v_segment, j_segment, scope, top_n)`** -
  Match single TCR
- **`match_tcr_many_df(db, cdr3, v, j, scope, top_n, progress)`** -
  Batch match with progress bar

**Search Scope Format:** `"s,i,d,t"` where: - `s` = max substitutions -
`i` = max insertions - `d` = max deletions - `t` = max total edit
distance

Examples: - `"0,0,0,0"` - Exact match only - `"2,1,2,3"` - Max 2
substitutions, 1 insertion, 2 deletions, 3 total edits - `"3,1,2,3"` -
More permissive fuzzy matching

### Seurat Integration

- **`vdj_attach_10x_vdj_v2(obj, vdj_dir, ...)`** - Attach 10x VDJ data
  to Seurat object
- **`vdj_attach_10x_vdj_v2_batch(obj, matchdf, vdj_dir_template)`** -
  Batch attach multiple samples
- **`vdj_collapse_pairs_seurat(obj, ...)`** - Collapse TRA/TRB chains
  into pairs

## Vignettes

Comprehensive guides covering all aspects of vdjmatchR:

### 📚 Getting Started

[**Getting started with
vdjmatchR**](https://furlan-lab.github.io/vdjmatchR/articles/getting-started.md)
Installation, basic usage, and your first TCR match

### 🚀 Batch Matching

[**Batch matching
clonotypes**](https://furlan-lab.github.io/vdjmatchR/articles/batch-matching.md)
High-throughput matching with parallel processing and progress bars

### 🎯 Scoring and Search

[**Scoring and search
scope**](https://furlan-lab.github.io/vdjmatchR/articles/scoring-and-scope.md)
Understanding fuzzy matching, edit distances, and scoring methods

### 💾 Database Management

[**Database management and
filtering**](https://furlan-lab.github.io/vdjmatchR/articles/database-management.md)
Loading, filtering, updating databases, and converting to data.table

### 📊 Analysis & Visualization

[**TCR motif
analysis**](https://furlan-lab.github.io/vdjmatchR/articles/tcr-motif-analysis.md)
Create sequence motif logos and visualize conserved patterns in
epitope-specific TCRs

### 🧬 Seurat Integration

[**Seurat v5 + 10x VDJ v2
Integration**](https://furlan-lab.github.io/vdjmatchR/articles/seurat-vdj-integration.md)
Complete workflow for integrating 10x VDJ data with Seurat single-cell
analysis

## Performance Features

### Parallel Processing

- Automatic multi-core processing via Rayon (Rust)
- Linear speedup with CPU cores (e.g., 8x faster on 8-core machines)
- No configuration required

### Progress Bars

- Real-time progress tracking for large datasets
- Automatic chunking for memory efficiency
- Summary statistics on completion

Example with 10,000 TCRs:

    Matching 10000 TCRs against database (parallel processing enabled)
    |=============================================| 100%
    Matched 10000 queries, found 3842 hits

### Memory Efficient

- Streaming database access
- Minimal memory overhead
- Scales to millions of queries

## Example Workflows

### Explore Database

``` r
library(vdjmatchR)
library(data.table)

# Load and convert database
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))
dt <- db_to_table(db)

# Summary
db_summary(db)

# Find most common epitopes
dt[, .N, by = antigen_epitope][order(-N)][1:10]

# Viral epitopes
dt[antigen_species %in% c("CMV", "EBV", "InfluenzaA"),
   .N, by = .(antigen_species, antigen_epitope)][order(-N)]
```

### Match Seurat TCRs

``` r
library(vdjmatchR)
library(Seurat)

# Load Seurat object and attach VDJ data
obj <- readRDS("seurat_object.rds")
obj <- vdj_attach_10x_vdj_v2(obj, vdj_dir = "path/to/vdj")

# Collapse to paired chains
obj <- vdj_collapse_pairs_seurat(obj, write_meta = TRUE)

# Load and filter database
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))
db_trb <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)

# Extract TCR sequences
md <- obj@meta.data
queries <- md[!is.na(md$vdj.trb.cdr3_aa), ]

# Match against VDJdb
hits <- match_tcr_many_df(
  db_trb,
  queries$vdj.trb.cdr3_aa,
  queries$vdj.trb.v_gene,
  queries$vdj.trb.j_gene,
  scope = "3,1,2,3",
  top_n = 1,
  progress = TRUE
)

# Add epitope annotations back to Seurat
obj$epitope <- NA
obj$epitope[match(hits$cell, rownames(obj@meta.data))] <- hits$epitope
```

## Documentation

- **Website:** <https://furlan-lab.github.io/vdjmatchR/>
- **Reference:** Complete function documentation with examples
- **Vignettes:** Step-by-step guides for common workflows

## Database Information

vdjmatchR uses the [VDJdb
database](https://github.com/antigenomics/vdjdb-db) - a curated
repository of TCR sequences with known antigen specificities.

**Bundled Databases:** - **Slim** (~2.3 MB compressed): Essential
columns only - **Fat** (~5.9 MB compressed): Complete metadata

**Database Contents:** - Human and mouse TCR sequences - Known epitope
specificities - V/J gene annotations - MHC restriction information -
Confidence scores - Publication references

## Contributing

Contributions welcome! Please:

1.  Fork the repository
2.  Create a feature branch
3.  Make your changes
4.  Add tests if applicable
5.  Submit a pull request

## Citation

If you use vdjmatchR in your research, please cite:

    vdjmatchR: Fast TCR sequence matching with VDJdb in R
    https://github.com/furlan-lab/vdjmatchR

And the VDJdb database:

    Shugay M et al. VDJdb: a curated database of T-cell receptor sequences
    with known antigen specificity. Nucleic Acids Research (2018).

## License

\[Add your license here\]

## Credits

- Built with [extendr](https://extendr.github.io/)
- Powered by [VDJdb](https://github.com/antigenomics/vdjdb-db)
- Parallel processing via [Rayon](https://github.com/rayon-rs/rayon)

------------------------------------------------------------------------

Made with ❤️ by the Furlan Lab
