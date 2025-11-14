# Getting started with vdjmatchR

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

## Quick tour

``` r
library(vdjmatchR)

# Use the packaged slim VDJdb bundled with the package
path <- vdjdb_packaged_path(use_fat_db = FALSE)
db <- vdjdb_open_file(path)

# Optionally filter to human TRB with min score >= 2
fdb <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)

# Match a single clonotype, return a data.frame
df <- match_tcr_df(
  fdb,
  cdr3 = "CASSLGQAYEQYF",
  v_segment = "TRBV12-3",
  j_segment = "TRBJ2-7",
  scope = "0,0,0,0",
  top_n = 5
)
print(df)

# Optionally, point to your own VDJdb file for this session
# vdjdb_set_user_db("/path/to/vdjdb.slim.txt.gz", use_fat_db = FALSE)
# db <- vdjdb_open(use_fat_db = FALSE)

# Batch matching of multiple clonotypes
cdr3s <- c("CASSLGQAYEQYF", "CASSIRSSYEQYF")
vs    <- c("TRBV12-3",      "TRBV7-9")
js    <- c("TRBJ2-7",       "TRBJ2-7")
res_many <- match_tcr_many_df(fdb, cdr3s, vs, js, scope = "0,0,0,0", top_n = 3)
head(res_many)
```

## Search scope

The `scope` argument controls fuzzy matching tolerance using the format
`"s,i,d,t"` (substitutions, insertions, deletions, total) or `"s,id,t"`
where insertions and deletions share the same limit. Use `"0,0,0,0"` for
exact matching.

``` r
# Allow up to 2 substitutions and 1 insertion/deletion, total up to 3 edits
scope <- "2,1,2,3"
```

## Performance tips

- Pre-filter the database (species/gene/score) to speed up matching.
- Use `top_n` to cap returned hits per query.
- For large batches, prefer
  [`match_tcr_many_df()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many_df.md).

## Reproducibility

Document the database version and filters used in your analysis. Store
the VDJdb file path used
([`vdjdb_path()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_path.md))
in your project metadata.
