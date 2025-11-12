# vdjmatchR Implementation Summary

Date: 2025-11-11

This document summarizes the work done to create the `vdjmatchR` repository, an R
package exposing the Rust TCR matching functionality (derived from `vdjmatch-rs`)
via `extendr`.

## Overview
- Scaffolded a new R package in `~/develop/vdjmatchR`.
- Embedded a Rust crate (`src/rust`) with core matching logic and `extendr` wrappers.
- Added R convenience functions returning data.frames for ergonomic usage.
- Implemented pkgdown documentation site configuration and multiple vignettes.
- Included a GitHub Actions workflow to build and deploy the pkgdown site.

## Repository Structure (key files)
- R package
  - `DESCRIPTION` — package metadata (updated for vignettes/pkgdown).
  - `NAMESPACE` — loads dynamic library; exports default patterns.
  - `R/zzz.R` — `.onLoad` to load shared library.
  - `R/match.R` — high-level R helpers (data.frame results, filters, downloads).
  - `R/RDatabase.R` — roxygen topic for the `RDatabase` class.
  - `README.md` — quickstart and examples.
  - `.Rbuildignore`, `.gitignore` — ignores incl. pkgdown artifacts.
  - `LICENSE`, `LICENSE.md` — MIT license metadata.
- Rust crate (extendr + core)
  - `src/rust/Cargo.toml` — cdylib crate with dependencies (`extendr-api`, etc.).
  - `src/rust/build.rs` — enables rextendr wrapper generation.
  - `src/rust/src/lib.rs` — extendr module and exports.
  - `src/rust/src/{alignment,database,error,filtering,matching,scoring,sequence,utils}.rs` — copied from `vdjmatch-rs`.
- R build glue
  - `src/Makevars`, `src/Makevars.win` — link flags via `rextendr`.
  - `src/entrypoint.c` — minimal C stub (no registration; extendr handles it).
- Pkgdown and docs
  - `_pkgdown.yml` — site configuration, navbar, and reference groups.
  - `.github/workflows/pkgdown.yaml` — CI to build and deploy site.
  - `vignettes/getting-started.Rmd` — intro vignette.
  - `vignettes/scoring-and-scope.Rmd` — scoring and tolerance.
  - `vignettes/batch-matching.Rmd` — batch matching workflows.
  - `vignettes/database-management.Rmd` — downloading and filtering DB.

## Rust API exposed to R (extendr)
File: `src/rust/src/lib.rs`

- Type
  - `RDatabase` (wrapper holding `Database`):
    - `new_from_file(path: &str) -> RDatabase`
    - `new_from_vdjdb(use_fat_db: bool) -> RDatabase` (ensures local VDJdb then loads)
    - `len(&self) -> i32`
    - `filter(&self, species: Option<String>, gene: Option<String>, min_vdjdb_score: i32) -> RDatabase`
    - `filter_by_epitope_size(&self, min_size: i32) -> RDatabase`
- Functions
  - `match_tcr(db, cdr3, v_segment, j_segment, scope, top_n) -> List` (single clonotype)
  - `match_tcr_many(db, cdr3[], v_segment[], j_segment[], scope, top_n) -> List` (batch)
  - `vdjdb_ensure(use_fat_db: bool) -> String` (returns TSV path; downloads if missing)
  - `vdjdb_update() -> ()` (downloads/refreshes slim and full TSVs)

Notes
- Scoring currently uses the simple mismatch approach defined in the copied `matching/scoring` modules.
- `scope` parsing supports `"s,i,d,t"` or `"s,id,t"` formats; `"0,0,0,0"` means exact.
- Batch results include `query_index` (1-based for R) and query metadata columns.

## R Convenience API
File: `R/match.R`

- `match_tcr_df(db, cdr3, v_segment = "", j_segment = "", scope = "0,0,0,0", top_n = 0)` → data.frame
- `match_tcr_many_df(db, cdr3[], v_segment[], j_segment[], scope = "0,0,0,0", top_n = 0)` → data.frame
- `filter_db(db, species = NULL, gene = NULL, min_vdjdb_score = 0)` → `RDatabase`
- `filter_db_by_epitope_size(db, min_size = 1)` → `RDatabase`
- `vdjdb_download(use_fat_db = FALSE)` → TSV file path
- `vdjdb_update_all()` → downloads/refreshes VDJdb

## Build and Installation
Requirements
- Rust ≥ 1.70 and cargo
- R toolchain with `rextendr` (for compile/link flags)

Install steps
- In R: `install.packages("rextendr")`
- In shell from package root: `R CMD INSTALL .`

## Pkgdown Site
- Local build: in R, run `roxygen2::roxygenise(); pkgdown::build_site()`.
- CI deploy: GitHub Actions builds `docs/` and publishes to Pages (enable Pages for the repo).

## Vignettes
- Getting started: walkthrough of install, download, filtering, and matching.
- Scoring and scope: guidance on `scope` and interpreting scores.
- Batch matching: efficient multi-query usage.
- Database management: ensuring/updating VDJdb and filtering strategies.

## Notable Implementation Details
- `src/entrypoint.c` contains only a dummy symbol to avoid conflicting with
  extendr’s auto-registration; the Rust side provides the `R_init_*` hook.
- Database download uses `reqwest` with blocking + rustls TLS; requires network
  when calling `vdjdb_ensure()` or `vdjdb_update()`.
- Core Rust modules are currently copied from `vdjmatch-rs` to keep this repo
  self-contained.

## Next Steps (optional)
- Replace placeholder URLs in `DESCRIPTION` and `_pkgdown.yml` with the real GitHub org/repo.
- Refactor to depend on a library crate version of `vdjmatch-rs` to avoid code duplication (workspace or git dependency).
- Add tests (R and Rust) and CI for build checks.
- Add more R-friendly result shaping (e.g., tibble, factor levels for segments).
- Expose additional configuration knobs (e.g., scoring mode toggle) through the R API.

