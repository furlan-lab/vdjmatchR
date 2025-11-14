# Attach 10x 5' VDJ v2 data to a Seurat v5 object

Reads Cell Ranger VDJ v2 outputs (filtered_contig_annotations.csv and
clonotypes.csv), flexibly maps 10x barcodes to Seurat cell names, and
stores heavy tables under `obj@tools[[tool_key]]`. Optionally writes
minimal metadata (10x clonotype id/size).

## Usage

``` r
vdj_attach_10x_vdj_v2(
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
)
```

## Arguments

- seurat:

  A Seurat v5 object.

- vdj_dir:

  Path to the 10x VDJ directory containing
  `filtered_contig_annotations.csv` and `clonotypes.csv`.

- sample_id:

  Optional sample identifier. If NULL, inferred from `vdj_dir` (parent
  folder name if `vdj_dir` ends with `vdj_t`).

- barcode_prefix:

  Optional prefix to help exact matching (tries `prefix + barcode` and
  `prefix + "_" + barcode`).

- cell_regex:

  Optional regex with one capture group to extract barcodes from cell
  names.

- cell_map:

  Optional explicit mapping (named character vector barcode-\>cell or a
  data.frame with `barcode` and `cell`). Overrides other strategies.

- tool_key:

  Name under `obj@tools` where results are stored (default `"vdj"`).

- add_meta:

  If TRUE, write `vdj.clone_id_10x` and `vdj.clone_size_10x` to metadata
  for mapped cells.

- include_cdr12:

  If TRUE, attempt to import CDR1 and CDR2 amino-acid sequences when
  columns `cdr1` and `cdr2` exist in the contigs file.

- include_nt_subseq:

  If TRUE, and when corresponding columns exist (e.g., `cdr3_nt`,
  `cdr1_nt`, `cdr2_nt`), also write nucleotide sequences for these
  subsequences.

- expected_t_fraction:

  Expected fraction of T cells among queried Seurat cells. Used to gate
  warnings about unmatched cells. A warning is emitted only if
  `(unmatched / queried) > (1 - expected_t_fraction)`. Default: `0.5`.

## Value

Updated Seurat object (invisible).

## Details

Assumptions: human, TRA/TRB only. No code path depends on Seurat being
attached; the function only uses base slots and data.frame operations to
keep Seurat as an optional dependency.

Mapping strategies (applied in priority order, filling only unmatched
barcodes):

- `cell_map`: explicit mapping (named character vector barcode-\>cell or
  data.frame with columns `barcode` and `cell`).

- `cell_regex`: regex with a single capture group applied to Seurat cell
  names to extract the barcode; unique matches are used.

- Exact match: Seurat cell name exactly equals the 10x barcode (e.g.,
  "AAAC...-1").

- `barcode_prefix`: attempt exact match using `paste0(prefix, barcode)`
  and `paste0(prefix, "_", barcode)`.

- Suffix match: Seurat cell name ends with the 10x barcode (common when
  names are prefixed).

Multi-sample: All stored tables include `sample_id`. Mapping does not
restrict to any metadata field; it operates purely on cell barcodes
(with optional prefixes/regex).
