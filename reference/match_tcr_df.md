# Match one clonotype and return a data.frame

Match one clonotype and return a data.frame

## Usage

``` r
match_tcr_df(
  db,
  cdr3,
  v_segment = "",
  j_segment = "",
  scope = "0,0,0,0",
  top_n = 0L
)
```

## Arguments

- db:

  an RDatabase object

- cdr3:

  CDR3 amino-acid sequence

- v_segment:

  V segment (optional; empty string to ignore)

- j_segment:

  J segment (optional; empty string to ignore)

- scope:

  search scope string like "0,0,0,0" or "2,1,2,3"

- top_n:

  keep top N hits (per query)

## Value

data.frame with matching hits
