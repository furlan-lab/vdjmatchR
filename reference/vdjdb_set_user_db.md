# Set a user-specified VDJdb path to use by default

Stores the path in R options for the session. Call with
`use_fat_db = TRUE` to set the full database, or FALSE for the slim one.

## Usage

``` r
vdjdb_set_user_db(path, use_fat_db = FALSE)
```

## Arguments

- path:

  Path to a VDJdb TSV/TSV.GZ file

- use_fat_db:

  logical; TRUE for full db, FALSE for slim (default)
