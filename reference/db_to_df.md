# Convert VDJdb database to data.frame

Extracts all entries from the database and returns as a data.frame. For
larger databases, consider using
[`db_to_table()`](https://furlan-lab.github.io/vdjmatchR/reference/db_to_table.md)
which returns a data.table for better performance.

## Usage

``` r
db_to_df(db)
```

## Arguments

- db:

  an RDatabase object

## Value

data.frame with database entries

## Examples

``` r
if (FALSE) { # \dontrun{
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = FALSE))
df <- db_to_df(db)
head(df)
} # }
```
