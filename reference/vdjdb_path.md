# Preferred VDJdb path (user-specified if set, else packaged)

Returns a user-specified path if configured via
[vdjdb_set_user_db](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_set_user_db.md),
otherwise returns the packaged file path.

## Usage

``` r
vdjdb_path(use_fat_db = FALSE)
```

## Arguments

- use_fat_db:

  logical; TRUE for full db, FALSE for slim (default)

## Value

string file path
