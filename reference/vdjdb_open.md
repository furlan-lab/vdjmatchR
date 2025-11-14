# Open the preferred VDJdb (user path if set, else packaged)

Open the preferred VDJdb (user path if set, else packaged)

## Usage

``` r
vdjdb_open(use_fat_db = FALSE)
```

## Arguments

- use_fat_db:

  logical; TRUE for full db, FALSE for slim (default)

## Value

An in-memory VDJdb handle created by the Rust backend.
