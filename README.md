# vdjmatchR

<p align="center">
  <img src="blob/logo2.png" alt="vdjmatchR logo" width="300">
</p>

R interface to a Rust implementation of TCR sequence matching (derived from `vdjmatch-rs`), built via [`extendr`](https://extendr.github.io/).

## Overview

Fast TCR sequence matching against the VDJdb database with a Rust-powered core.

### Key Features

✨ **Fast Matching** - Parallel processing via Rayon for high throughput
🔍 **Fuzzy Matching** - Configurable edit distance tolerance
📊 **Database Tools** - Easy filtering and conversion to data.table
🧬 **Seurat Integration** - Seamless integration with 10x VDJ data
📈 **Progress Tracking** - Real-time progress bars for large datasets

## Documentation

📖 **Website:** https://furlan-lab.github.io/vdjmatchR/
📚 **Vignettes:** Complete guides for all workflows
💡 **Examples:** See [index.md](index.md) for detailed usage examples

## Quickstart
1. Ensure Rust (>= 1.70) and an R toolchain are installed.
2. In R, install `rextendr` to aid linking (optional but recommended):
   `install.packages("rextendr")`
3. Build + install the package from the package root:
   `R CMD INSTALL .`

## Examples
```r
library(vdjmatchR)

# Use the packaged slim VDJdb (bundled with the package)
path <- vdjdb_packaged_path(use_fat_db = FALSE)
db <- vdjdb_open_file(path)

# Optionally update the packaged database into the package's extdata and use it
# updated_path <- vdjdb_update_latest(use_fat_db = FALSE)
# db <- vdjdb_open_file(updated_path)

# Or point to your own database file for this session
# vdjdb_set_user_db("/path/to/vdjdb.slim.txt.gz", use_fat_db = FALSE)
# path <- vdjdb_path(FALSE)
# db <- vdjdb_open_file(path)

# Or use the packaged full ("fat") VDJdb
# path_full <- vdjdb_packaged_path(use_fat_db = TRUE)
# db <- vdjdb_open_file(path_full)

# Filter for human TRB and higher-confidence entries
fdb <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)

# Match a single clonotype; get a data.frame
df <- match_tcr_df(fdb, cdr3 = "CASSLGQAYEQYF", v_segment = "TRBV12-3", j_segment = "TRBJ2-7", scope = "0,0,0,0", top_n = 5)
print(df)

# Batch match: vectors of clonotypes
cdr3s <- c("CASSLGQAYEQYF", "CASSIRSSYEQYF")
vs    <- c("TRBV12-3", "TRBV7-9")
js    <- c("TRBJ2-7",  "TRBJ2-7")
res_many <- match_tcr_many_df(fdb, cdr3s, vs, js, scope = "0,0,0,0", top_n = 3)
head(res_many)
```

## R API
- `vdjdb_open_file(path)` / `vdjdb_open_packaged(use_fat_db = FALSE)` / `vdjdb_open(use_fat_db = FALSE)`
- `vdjdb_len(db)`
- `filter_db(db, species = NULL, gene = NULL, min_vdjdb_score = 0)`
- `filter_db_by_epitope_size(db, min_size = 1)`
- `match_tcr(db, cdr3, v_segment, j_segment, scope, top_n)` → list
- `match_tcr_df(...)` → data.frame
- `match_tcr_many(db, cdr3, v_segment, j_segment, scope, top_n)` → list
- `match_tcr_many_df(...)` → data.frame
- `vdjdb_packaged_path(use_fat_db = FALSE)` / `vdjdb_path(use_fat_db = FALSE)`
- `vdjdb_set_user_db(path, use_fat_db = FALSE)` (configure your own DB file)

## Vignettes

Comprehensive guides covering all aspects of vdjmatchR:

- 📚 **[Getting started](https://furlan-lab.github.io/vdjmatchR/articles/getting-started.html)** - Installation, basic usage, first TCR match
- 🚀 **[Batch matching](https://furlan-lab.github.io/vdjmatchR/articles/batch-matching.html)** - High-throughput matching with parallel processing
- 🎯 **[Scoring and search](https://furlan-lab.github.io/vdjmatchR/articles/scoring-and-scope.html)** - Understanding fuzzy matching and scoring
- 💾 **[Database management](https://furlan-lab.github.io/vdjmatchR/articles/database-management.html)** - Loading, filtering, and exploring databases
- 🧬 **[Seurat integration](https://furlan-lab.github.io/vdjmatchR/articles/seurat-vdj-integration.html)** - Complete 10x VDJ workflow

## Notes

- Results are returned as data.frames via R wrappers for convenience; core Rust functions return plain lists for flexibility.
- Use `scope` like `"0,0,0,0"` for exact or e.g. `"2,1,2,3"` for controlled fuzzy matching.
- Parallel processing is enabled automatically for batch matching.
- Progress bars appear for large datasets (>5000 queries).
