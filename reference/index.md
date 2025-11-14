# Package index

## TCR Matching

Core functions for matching TCR sequences against VDJdb

- [`match_tcr()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr.md)
  : Match a single clonotype against the database. Returns a list of
  columns (vector-of-equal-length) suitable for as.data.frame in R.
- [`match_tcr_df()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_df.md)
  : Match one clonotype and return a data.frame
- [`match_tcr_many()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many.md)
  : Batch match: vectors of cdr3/v/j; returns stacked results with query
  metadata. Uses parallel processing via Rayon for improved performance.
- [`match_tcr_many_df()`](https://furlan-lab.github.io/vdjmatchR/reference/match_tcr_many_df.md)
  : Match many clonotypes and return a data.frame stacked across queries

## Database Loading & Management

Functions to load, filter, and manage VDJdb databases

- [`vdjdb_open_file()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_open_file.md)
  : Open a VDJdb TSV/TSV.GZ via the Rust backend.
- [`vdjdb_open()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_open.md)
  : Open the preferred VDJdb (user path if set, else packaged)
- [`vdjdb_open_packaged()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_open_packaged.md)
  : Open the packaged VDJdb via the Rust backend
- [`vdjdb_packaged_path()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_packaged_path.md)
  : Locate packaged VDJdb file bundled with vdjmatchR
- [`vdjdb_path()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_path.md)
  : Preferred VDJdb path (user-specified if set, else packaged)
- [`vdjdb_set_user_db()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_set_user_db.md)
  : Set a user-specified VDJdb path to use by default
- [`vdjdb_len()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_len.md)
  : Number of rows stored in the in-memory VDJdb handle.
- [`vdjdb_download()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_download.md)
  : Back-compat alias for updating/downloading VDJdb (one variant)
- [`vdjdb_ensure()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_ensure.md)
  : Ensure VDJdb exists locally and return the path.
- [`vdjdb_update()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update.md)
  : Download/update the VDJdb files (slim and fat).
- [`vdjdb_update_all()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update_all.md)
  : Update/download both slim and fat VDJdb files (disabled)
- [`vdjdb_update_latest()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update_latest.md)
  : Update (download) the VDJdb file and return its path into extdata
- [`vdjdb_ensure_into()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_ensure_into.md)
  : Ensure VDJdb exists in the specified directory and return the path.
- [`vdjdb_update_into()`](https://furlan-lab.github.io/vdjmatchR/reference/vdjdb_update_into.md)
  : Download/update the VDJdb files (slim and fat) into the specified
  directory.

## Database Filtering

Filter databases by species, gene, score, and epitope representation

- [`filter_db()`](https://furlan-lab.github.io/vdjmatchR/reference/filter_db.md)
  : Filter database entries by species, gene, and minimum VDJdb score.
- [`filter_db_by_epitope_size()`](https://furlan-lab.github.io/vdjmatchR/reference/filter_db_by_epitope_size.md)
  : Filter by minimum epitope size (unique CDR3 per epitope).

## Database Conversion & Exploration

Convert databases to R data structures for easy manipulation

- [`db_to_table()`](https://furlan-lab.github.io/vdjmatchR/reference/db_to_table.md)
  : Convert VDJdb database to data.table
- [`db_to_df()`](https://furlan-lab.github.io/vdjmatchR/reference/db_to_df.md)
  : Convert VDJdb database to data.frame
- [`db_summary()`](https://furlan-lab.github.io/vdjmatchR/reference/db_summary.md)
  : Quick database summary

## Seurat Integration

Functions for integrating 10x VDJ data with Seurat objects

- [`vdj_attach_10x_vdj_v2()`](https://furlan-lab.github.io/vdjmatchR/reference/vdj_attach_10x_vdj_v2.md)
  : Attach 10x 5' VDJ v2 data to a Seurat v5 object
- [`vdj_attach_10x_vdj_v2_batch()`](https://furlan-lab.github.io/vdjmatchR/reference/vdj_attach_10x_vdj_v2_batch.md)
  : Attach multiple 10x 5' VDJ v2 runs using a caps-token template and
  prefix map
- [`vdj_collapse_pairs_seurat()`](https://furlan-lab.github.io/vdjmatchR/reference/vdj_collapse_pairs_seurat.md)
  : Collapse per-cell TRA/TRB into 1 TRB + up to 2 TRA pairs

## Package Information

- [`vdjmatchR-package`](https://furlan-lab.github.io/vdjmatchR/reference/vdjmatchR-package.md)
  : vdjmatchR: R Interface to VDJmatch (Rust)
