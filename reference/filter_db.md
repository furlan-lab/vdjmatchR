# Filter database via species/gene/score and return a new RDatabase

Filter database via species/gene/score and return a new RDatabase

## Usage

``` r
filter_db(db, species = NULL, gene = NULL, min_vdjdb_score = 0L)
```

## Arguments

- db:

  an RDatabase object

- species:

  optional species filter (e.g. "HomoSapiens")

- gene:

  optional gene filter (e.g. "TRB")

- min_vdjdb_score:

  minimum VDJdb score (0-3)

## Value

an RDatabase object
