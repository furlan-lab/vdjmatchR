# Batch matching clonotypes

``` r
#| echo: false
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", eval = FALSE)
```

This vignette demonstrates how to match multiple clonotypes efficiently
using
[`match_tcr_many_df()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many_df.md).

## Setup

``` r
library(vdjmatchR)
path <- vdjdb_packaged_path(FALSE)
db <- RDatabase$new_from_file(path)
fdb <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)
```

## Prepare inputs

``` r
cdr3s <- c("CASSLGQAYEQYF", "CASSIRSSYEQYF", "CASSPGQGYEQYF")
vs    <- c("TRBV12-3",      "TRBV7-9",       "TRBV20-1")
js    <- c("TRBJ2-7",       "TRBJ2-7",       "TRBJ2-5")
```

## Run batch matching

``` r
res <- match_tcr_many_df(fdb, cdr3s, vs, js, scope = "0,0,0,0", top_n = 3)
head(res)
```

Each result row includes the 1-based `query_index` that maps back to the
input vectors, along with the query clonotype and hit columns.

## Tips

- Ensure the three input vectors have equal length.
- Use `top_n` to cap returned hits per query to a manageable size.
- Pre-filter the database to reduce compute and noise.
