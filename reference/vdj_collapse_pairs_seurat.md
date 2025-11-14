# Collapse per-cell TRA/TRB into 1 TRB + up to 2 TRA pairs

Selects one TRB per cell and up to two TRA chains ordered by support,
then emits a `pairs` table in `obj@tools[[tool_key]]$pairs` with one or
two rows per cell: the dominant pair (TRA1+TRB) and, if present and
requested, the secondary pair (TRA2+TRB).

## Usage

``` r
vdj_collapse_pairs_seurat(
  seurat,
  tool_key = "vdj",
  drop_multi_trb = TRUE,
  keep_two_alpha = TRUE,
  alpha_support = c("umis", "reads"),
  min_umis = 0L,
  min_reads = 0L,
  require_productive = TRUE,
  require_full_length = FALSE,
  require_high_confidence = FALSE,
  write_meta = TRUE,
  update_primary_meta = TRUE
)
```

## Arguments

- seurat:

  Seurat v5 object previously processed by
  `vdj_attach_10x_vdj_v2[_batch]`.

- tool_key:

  Character; tools slot name (default "vdj").

- drop_multi_trb:

  Logical; drop cells with \>1 TRB after QC (default TRUE). If FALSE,
  picks the best TRB.

- keep_two_alpha:

  Logical; if TRUE, produce a secondary pair using TRA2+TRB when
  available (default TRUE).

- alpha_support:

  Character; support metric for ranking alphas: "umis" (default) or
  "reads".

- min_umis:

  Integer; discard contigs below this UMI count (default 0).

- min_reads:

  Integer; discard contigs below this read count (default 0).

- require_productive:

  Logical; require productive contigs (default TRUE).

- require_full_length:

  Logical; require full_length when column exists (default FALSE).

- require_high_confidence:

  Logical; require high_confidence when column exists (default FALSE).

- write_meta:

  Logical; write metadata columns described above (default TRUE).

- update_primary_meta:

  Logical; if TRUE, also set `vdj.tra.*` and `vdj.trb.*` to the primary
  pair (default TRUE).

## Value

Updated Seurat object (invisible).

## Details

Selection heuristic:

- Filter contigs by QC: `productive`, `full_length`, `high_confidence`
  (configurable).

- TRB: select exactly one best contig by UMIs, then reads, then
  high_confidence. If more than one TRB exists and
  `drop_multi_trb = TRUE`, the cell is dropped.

- TRA: order candidate alphas by UMIs (or reads) descending, then reads,
  then high_confidence; keep the top one (dominant) and optionally the
  second.

Metadata:

- Writes `vdj.tra1.*`, `vdj.tra2.*` (if present), and `vdj.trb.*` fields
  per cell when `write_meta = TRUE`. Also writes flags `vdj.has_tra2`,
  `vdj.multi_trb`, and `vdj.pair_qc`.

- Optionally updates the primary `vdj.tra.*` and `vdj.trb.*` fields to
  TRA1/TRB.
