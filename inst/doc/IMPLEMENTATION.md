# vdjmatchR — Implementation Notes (Current State)

This document describes the current implementation and public surface of the
`vdjmatchR` R package, including TCR matching against VDJdb and Seurat + 10x 5' VDJ v2 integration.

## Overview
- Exposes TCR matching via Rust/extendr bindings (no network I/O).
- Provides Seurat utilities to ingest Cell Ranger 10x 5' VDJ v2 outputs into a Seurat v5 object.
- Stores heavy VDJ tables under `obj@tools$vdj` and writes concise, per-cell VDJ metadata.

## Repository Structure (key files)
- R API
  - `R/match.R` — high-level matching helpers (data.frame results, filters).
  - `R/vdjdb.R` — packaged DB paths and user DB configuration.
  - `R/seurat_vdj.R` — Seurat + 10x VDJ integration functions and diagnostics.
- Rust + build glue
  - `src/` — extendr integration and static library glue (platform-specific Makevars).
- Docs
  - `vignettes/` — demonstration and integration workflows.
  - `_pkgdown.yml` — site configuration.

## Matching API (R)
- `match_tcr_df(db, cdr3, v_segment = "", j_segment = "", scope = "0,0,0,0", top_n = 0)` → data.frame
- `match_tcr_many_df(db, cdr3[], v_segment[], j_segment[], scope = "0,0,0,0", top_n = 0)` → data.frame
- Database helpers: `filter_db`, `filter_db_by_epitope_size`, `vdjdb_packaged_path`,
  `vdjdb_path`, `vdjdb_set_user_db`, `vdjdb_update_latest` (no network), `vdjdb_download` (no network).

## Seurat + 10x VDJ Integration (R/seurat_vdj.R)

### Single-run attach
`vdj_attach_10x_vdj_v2(
  seurat,
  vdj_dir,
  sample_id = NULL,
  barcode_prefix = NULL,
  cell_regex = NULL,
  cell_map = NULL,
  tool_key = "vdj",
  add_meta = TRUE,
  include_cdr12 = FALSE,
  include_nt_subseq = FALSE,
  expected_t_fraction = 0.5
)`

- Input: Cell Ranger VDJ v2 outputs in `vdj_dir` (`filtered_contig_annotations.csv`, `clonotypes.csv`).
- Sample: `sample_id` inferred from the path if missing.
- Mapping (barcode-based; no metadata gating):
  - Priority: explicit `cell_map` → `cell_regex` capture → exact equality → prefix-assisted equality (`prefix + barcode`, `prefix + "_" + barcode`) → suffix fallback.
  - Per-run mapping table: 10x `barcode` → Seurat `cell` with `strategy` (unique matches only).
- Heavy tables: append to `obj@tools[[tool_key]]$contigs` and `$clonotypes_raw` (adds `sample_id`).
- Per-cell minimal metadata:
  - `vdj.clone_id_10x`: majority of `raw_clonotype_id` across a cell’s contigs; tie-break by summed UMIs.
  - `vdj.clone_size_10x`: frequency from `clonotypes.csv` for the chosen clonotype (scoped by sample).
- Per-cell, per-chain metadata (best contig per (cell, TRA/TRB)):
  - Selection priority: productive==TRUE (if present) → highest UMIs → highest reads (`read_count` or `reads`) → high_confidence (if present).
  - Genes: `vdj.<tra|trb>.v_gene`, `vdj.<tra|trb>.j_gene`, `vdj.<tra|trb>.c_gene`.
  - Sequences: `vdj.<tra|trb>.cdr3_aa`; if `include_nt_subseq=TRUE` and present, `vdj.<tra|trb>.cdr3_nt`.
  - Counts: `vdj.<tra|trb>.umis`, `vdj.<tra|trb>.reads`.
  - Optional CDR1/2 AA: `vdj.<tra|trb>.cdr1_aa`, `vdj.<tra|trb>.cdr2_aa` when `include_cdr12=TRUE` and columns exist; `*_nt` when `include_nt_subseq=TRUE` and present.
- Diagnostics:
  - Queried cells: all cells if no prefix; else cells whose names start with `barcode_prefix`.
  - Unmatched fraction: (queried − matched) / queried.
  - Warning only if `unmatched_fraction > (1 - expected_t_fraction)`; default `expected_t_fraction=0.5`.
  - Summary row stored under `obj@tools[[tool_key]]$mapping_summary[[sample_id]]` with `queried`, `matched`, `unmatched`, `unmatched_fraction`.
- Tools structure also records `config$last_attach` (inputs and options).

### Batch attach with caps-token templates
`vdj_attach_10x_vdj_v2_batch(
  seurat,
  prefix_folders,
  vdj_dir_template,
  tool_key = "vdj",
  add_meta = TRUE,
  include_cdr12 = FALSE,
  include_nt_subseq = FALSE,
  expected_t_fraction = 0.5,
  placeholder = "FOLDERNAME",
  prefixes_col = c("prefixes","PREFIXES"),
  folders_col = "folders",
  cell_regex = NULL,
  cell_map = NULL
)`

- `prefix_folders`: data.frame with a prefix column (`PREFIXES` or `prefixes`) and any number of UPPERCASE columns referenced by the template (e.g., `FOLDERNAME1`, `FOLDERNAME2`).
- `vdj_dir_template`: path containing UPPERCASE tokens; replaced per-row with same-named columns. If no tokens are found, falls back to single `placeholder` substitution using `folders_col`.
- For each row: resolves `vdj_dir`, uses row’s prefix as `barcode_prefix`, and calls the single-run attach; results are appended.
- Missing directories are warned and skipped.

## Data written to Seurat
- `obj@tools[[tool_key]]`:
  - `$contigs`: rbind of per-run contigs with `sample_id`.
  - `$clonotypes_raw`: rbind of per-run clonotypes with `sample_id` and `clonotype_uid` (`sample:clonotype_id`).
  - `$mapping`: list of data.frames keyed by `sample_id` (10x `barcode`, mapped `cell`, `strategy`).
  - `$mapping_summary`: list of per-run data.frames with `queried`, `matched`, `unmatched`, `unmatched_fraction`.
  - `$config$last_attach`: inputs and options of the latest attach.
- `obj@meta.data`:
  - `vdj.clone_id_10x`, `vdj.clone_size_10x`.
  - `vdj.tra.*`, `vdj.trb.*` as described above.

## Vignettes
- `vignettes/seurat-vdj-integration.Rmd` demonstrates batch attach with caps-token templates, diagnostics via `mapping_summary`, and the per-chain metadata fields.
- Other vignettes: getting started, scoring/scope, batch matching, database management.

## Build and Installation
- Requires Rust toolchain and `rextendr` for linking; install via `R CMD INSTALL .` from the package root.

## Notes / Assumptions
- Assumes human TCR (TRA/TRB) in the Seurat helpers.
- Mapping is purely based on cell-barcode string operations and optional explicit/regex maps.
- No dependency on Seurat’s functions; uses base slots (`meta.data`, `tools`).

## Future Options
- Expose contig selection knobs (e.g., prioritize reads over UMIs, require `full_length`/`high_confidence`).
- Optional stricter mapping modes (e.g., disable suffix fallback).

