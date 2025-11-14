# Convert VDJdb database to data.table

Extracts all entries from the database and returns as a data.table for
easy inspection and manipulation in R.

## Usage

``` r
db_to_table(db)
```

## Arguments

- db:

  an RDatabase object

## Value

data.table with database entries

## Examples

``` r
if (FALSE) { # \dontrun{
# Load database
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))

# Convert to data.table
dt <- db_to_table(db)

# Filter human TRB entries
db_filtered <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 0)
dt_trb <- db_to_table(db_filtered)

# Now use data.table operations
library(data.table)
dt_trb[, .N, by = antigen_epitope][order(-N)][1:10]  # Top 10 epitopes
} # }
```
