# vdjmatchR

R interface to a Rust implementation of TCR sequence matching (derived from `vdjmatch-rs`), built via [`extendr`](https://extendr.github.io/).

## Overview
- Ships core matching logic in Rust for performance.
- Exposes an R API using extendr.
- Designed to work with VDJdb-formatted TSVs (slim or full).

## Quickstart
1. Ensure Rust (>= 1.70) and an R toolchain are installed.
2. In R, install `rextendr` to aid linking (optional but recommended):
   `install.packages("rextendr")`
3. Build + install the package from the package root:
   `R CMD INSTALL .`

## Examples
```r
library(vdjmatchR)

# Option A: Ensure/download VDJdb and load
path <- vdjdb_download(use_fat_db = FALSE)
db <- RDatabase$new_from_file(path)

# Option B: Load an existing TSV
# db <- RDatabase$new_from_file("/path/to/vdjdb.slim.txt")

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
- `RDatabase$new_from_file(path)` / `RDatabase$new_from_vdjdb(use_fat_db=FALSE)`
- `RDatabase$len()`
- `RDatabase$filter(species = NULL, gene = NULL, min_vdjdb_score = 0)`
- `RDatabase$filter_by_epitope_size(min_size = 1)`
- `match_tcr(db, cdr3, v_segment, j_segment, scope, top_n)` → list
- `match_tcr_df(...)` → data.frame
- `match_tcr_many(db, cdr3, v_segment, j_segment, scope, top_n)` → list
- `match_tcr_many_df(...)` → data.frame
- `vdjdb_download(use_fat_db = FALSE)` → file path
- `vdjdb_update_all()`

## Notes
- Results are returned as data.frames via R wrappers for convenience; core Rust functions return plain lists for flexibility.
- Use `scope` like `"0,0,0,0"` for exact or e.g. `"2,1,2,3"` for controlled fuzzy matching.

For a guided walkthrough, see the Getting started vignette: `vignettes/getting-started.Rmd`.
