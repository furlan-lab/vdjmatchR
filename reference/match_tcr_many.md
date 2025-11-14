# Batch match: vectors of cdr3/v/j; returns stacked results with query metadata. Uses parallel processing via Rayon for improved performance.

Batch match: vectors of cdr3/v/j; returns stacked results with query
metadata. Uses parallel processing via Rayon for improved performance.

## Usage

``` r
match_tcr_many(db, cdr3, v_segment, j_segment, scope, top_n)
```
