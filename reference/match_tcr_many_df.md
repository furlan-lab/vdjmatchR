# Match many clonotypes and return a data.frame stacked across queries

Match many clonotypes and return a data.frame stacked across queries

## Usage

``` r
match_tcr_many_df(
  db,
  cdr3,
  v_segment,
  j_segment,
  scope = "0,0,0,0",
  top_n = 0L
)
```

## Arguments

- db:

  an RDatabase object

- cdr3:

  character vector of CDR3 sequences

- v_segment:

  character vector of V segments (same length)

- j_segment:

  character vector of J segments (same length)

- scope:

  search scope string like "0,0,0,0" or "2,1,2,3"

- top_n:

  keep top N hits per query

## Value

data.frame with query metadata and hit columns
