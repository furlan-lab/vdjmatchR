# Attach multiple 10x 5' VDJ v2 runs using a caps-token template and prefix map

Convenience wrapper over
[`vdj_attach_10x_vdj_v2()`](https://furlan-lab.github.io/vdjmatchR/reference/vdj_attach_10x_vdj_v2.md)
to batch-attach many Cell Ranger VDJ outputs to a single Seurat object
when cell names were prefixed to avoid barcode collisions. The user
supplies a data.frame with a prefix column and any number of UPPERCASE
columns that correspond to tokens in a path template. All UPPERCASE
tokens present in the template (e.g., `FOLDERNAME1`, `FOLDERNAME2`) are
replaced per-row by the values of columns of the same names to locate
each run's VDJ output directory.

## Usage

``` r
vdj_attach_10x_vdj_v2_batch(
  seurat,
  prefix_folders,
  vdj_dir_template,
  tool_key = "vdj",
  add_meta = TRUE,
  include_cdr12 = FALSE,
  include_nt_subseq = FALSE,
  expected_t_fraction = 0.5,
  placeholder = "FOLDERNAME",
  prefixes_col = c("prefixes", "PREFIXES"),
  folders_col = "folders",
  cell_regex = NULL,
  cell_map = NULL
)
```

## Arguments

- seurat:

  A Seurat v5 object.

- prefix_folders:

  data.frame with at least a prefix column and one column per UPPERCASE
  token present in `vdj_dir_template`.

- vdj_dir_template:

  Character scalar path containing UPPERCASE token(s) such as
  `FOLDERNAME1`, `FOLDERNAME2` to be substituted from `prefix_folders`
  per row.

- tool_key:

  Name under `obj@tools` to store results (default `"vdj"`).

- add_meta:

  If TRUE, each run writes `vdj.clone_id_10x` and `vdj.clone_size_10x`
  for mapped cells.

- include_cdr12:

  If TRUE, forward to per-run attach to import CDR1/2 AA where
  available.

- include_nt_subseq:

  If TRUE, forward to per-run attach to include nucleotide sequences for
  CDR1/2/3 where available.

- expected_t_fraction:

  Expected fraction of T cells among queried Seurat cells; forwarded to
  per-run attach. A warning is emitted only if
  `(unmatched / queried) > (1 - expected_t_fraction)`. Default: `0.5`.

- placeholder:

  Token to replace in the template when no UPPERCASE tokens are detected
  (default `"FOLDERNAME"`).

- prefixes_col:

  Column name(s) in `prefix_folders` to look for the prefix. First found
  is used (default tries `"prefixes"`, then `"PREFIXES"`).

- folders_col:

  Column name in `prefix_folders` for folder names used only in
  single-placeholder fallback mode (default `"folders"`).

- cell_regex:

  Optional regex with one capture group to extract barcodes from cell
  names (forwarded to per-run attach).

- cell_map:

  Optional explicit mapping (forwarded to per-run attach).

## Value

Updated Seurat object (invisible).

## Details

Example template with two tokens:
`/path/to/MMCarT_240212_reseq/FOLDERNAME1/data/FOLDERNAME2/outs/per_sample_outs/FOLDERNAME2/vdj_t`

If no UPPERCASE tokens are detected in the template, the function falls
back to single-placeholder mode, replacing `placeholder` with values
from `folders_col`.
