# Package index

## Core functions

- [`match_tcr()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr.md)
  : Match a single clonotype against the database. Returns a list of
  columns (vector-of-equal-length) suitable for as.data.frame in R.
- [`match_tcr_df()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_df.md)
  : Match one clonotype and return a data.frame
- [`match_tcr_many()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many.md)
  : Batch match: vectors of cdr3/v/j; returns stacked results with query
  metadata.
- [`match_tcr_many_df()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many_df.md)
  : Match many clonotypes and return a data.frame stacked across queries
- [`vdjdb_download()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_download.md)
  : Back-compat alias for updating/downloading VDJdb (one variant)
- [`vdjdb_ensure()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_ensure.md)
  : Ensure VDJdb exists locally and return the path.
- [`vdjdb_packaged_path()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_packaged_path.md)
  : Locate packaged VDJdb file bundled with vdjmatchR
- [`vdjdb_path()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_path.md)
  : Preferred VDJdb path (user-specified if set, else packaged)
- [`vdjdb_set_user_db()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_set_user_db.md)
  : Set a user-specified VDJdb path to use by default
- [`vdjdb_update()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update.md)
  : Download/update the VDJdb files (slim and fat).
- [`vdjdb_update_all()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update_all.md)
  : Update/download both slim and fat VDJdb files (disabled)
- [`vdjdb_update_latest()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update_latest.md)
  : Update (download) the VDJdb file and return its path into extdata

## Database object

- [`RDatabase`](https://furlan-lab.github.io/vdjmatchR/reference/RDatabase.md)
  : RDatabase: Handle to a loaded VDJdb

## Helpers

- [`filter_db()`](https://furlan-lab.github.io/vdjmatchR/reference/filter_db.md)
  : Filter database via species/gene/score and return a new RDatabase
- [`filter_db_by_epitope_size()`](https://furlan-lab.github.io/vdjmatchR/reference/filter_db_by_epitope_size.md)
  : Filter database by minimum epitope size

## Package

- [`vdjmatchR-package`](https://furlan-lab.github.io/vdjmatchR/reference/vdjmatchR-package.md)
  : vdjmatchR: R Interface to VDJmatch (Rust)
