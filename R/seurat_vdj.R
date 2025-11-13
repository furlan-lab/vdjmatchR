#' Attach 10x 5' VDJ v2 data to a Seurat v5 object
#'
#' Reads Cell Ranger VDJ v2 outputs (filtered_contig_annotations.csv and clonotypes.csv),
#' flexibly maps 10x barcodes to Seurat cell names, and stores heavy tables under
#' `obj@tools[[tool_key]]`. Optionally writes minimal metadata (10x clonotype id/size).
#'
#' Assumptions: human, TRA/TRB only. No code path depends on Seurat being attached; the
#' function only uses base slots and data.frame operations to keep Seurat as an optional dependency.
#'
#' Mapping strategies (applied in priority order, filling only unmatched barcodes):
#' - `cell_map`: explicit mapping (named character vector barcode->cell or data.frame with
#'   columns `barcode` and `cell`).
#' - `cell_regex`: regex with a single capture group applied to Seurat cell names to extract
#'   the barcode; unique matches are used.
#' - Exact match: Seurat cell name exactly equals the 10x barcode (e.g., "AAAC...-1").
#' - `barcode_prefix`: attempt exact match using `paste0(prefix, barcode)` and `paste0(prefix, "_", barcode)`.
#' - Suffix match: Seurat cell name ends with the 10x barcode (common when names are prefixed).
#'
#' Multi-sample: All stored tables include `sample_id`. Mapping does not restrict to
#' any metadata field; it operates purely on cell barcodes (with optional prefixes/regex).
#'
#' @param seurat A Seurat v5 object.
#' @param vdj_dir Path to the 10x VDJ directory containing `filtered_contig_annotations.csv` and `clonotypes.csv`.
#' @param sample_id Optional sample identifier. If NULL, inferred from `vdj_dir` (parent folder name if `vdj_dir` ends with `vdj_t`).
#' @param barcode_prefix Optional prefix to help exact matching (tries `prefix + barcode` and `prefix + "_" + barcode`).
#' @param cell_regex Optional regex with one capture group to extract barcodes from cell names.
#' @param cell_map Optional explicit mapping (named character vector barcode->cell or a data.frame with `barcode` and `cell`). Overrides other strategies.
#' @param tool_key Name under `obj@tools` where results are stored (default `"vdj"`).
#' @param add_meta If TRUE, write `vdj.clone_id_10x` and `vdj.clone_size_10x` to metadata for mapped cells.
#' @param include_cdr12 If TRUE, attempt to import CDR1 and CDR2 amino-acid sequences when columns `cdr1` and `cdr2` exist in the contigs file.
#' @param include_nt_subseq If TRUE, and when corresponding columns exist (e.g., `cdr3_nt`, `cdr1_nt`, `cdr2_nt`), also write nucleotide sequences for these subsequences.
#' @param expected_t_fraction Expected fraction of T cells among queried Seurat cells. Used to gate warnings
#'   about unmatched cells. A warning is emitted only if `(unmatched / queried) > (1 - expected_t_fraction)`.
#'   Default: `0.5`.
#'
#' @return Updated Seurat object (invisible).
#' @export
vdj_attach_10x_vdj_v2 <- function(
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
) {
  # Basic checks without hard dependency on Seurat
  stopifnot(!is.null(seurat))
  if (!isS4(seurat) || !all(c("meta.data", "tools") %in% slotNames(seurat))) {
    stop("Expected a Seurat v5 object with slots 'meta.data' and 'tools'.")
  }

  # Resolve files
  vdj_dir <- normalizePath(vdj_dir, mustWork = TRUE)
  contigs_fp <- file.path(vdj_dir, "filtered_contig_annotations.csv")
  clonotypes_fp <- file.path(vdj_dir, "clonotypes.csv")
  if (!file.exists(contigs_fp)) stop(sprintf("Missing file: %s", contigs_fp))
  if (!file.exists(clonotypes_fp)) stop(sprintf("Missing file: %s", clonotypes_fp))

  # Infer sample_id if not provided
  if (is.null(sample_id) || !nzchar(sample_id)) {
    base <- basename(vdj_dir)
    if (base %in% c("vdj_t", "vdj_b", "vdj_t_gd")) {
      sample_id <- basename(dirname(vdj_dir))
    } else {
      sample_id <- base
    }
  }

  # Read 10x VDJ tables
  contigs <- utils::read.csv(contigs_fp, header = TRUE, stringsAsFactors = FALSE)
  clonotypes <- utils::read.csv(clonotypes_fp, header = TRUE, stringsAsFactors = FALSE)
  # Normalize expected columns exist
  required_cols <- c("barcode", "chain", "v_gene", "j_gene", "cdr3", "cdr3_nt")
  missing_cols <- setdiff(required_cols, names(contigs))
  if (length(missing_cols) > 0) {
    stop(sprintf("filtered_contig_annotations.csv is missing required columns: %s", paste(missing_cols, collapse = ", ")))
  }
  contigs$sample_id <- sample_id

  # Candidate cells: all cells in the object
  all_cells <- colnames(seurat)
  if (is.null(all_cells) || length(all_cells) == 0) all_cells <- rownames(seurat@meta.data)
  candidate_cells <- all_cells
  # If a barcode_prefix is provided (common in multi-sample objects), restrict
  # candidate cells to those belonging to this prefix. This prevents suffix
  # matching from mapping barcodes to cells in other samples that share the
  # same 10x barcode string.
  if (!is.null(barcode_prefix) && nzchar(barcode_prefix)) {
    candidate_cells <- candidate_cells[startsWith(candidate_cells, barcode_prefix)]
  }

  # Build mapping from 10x barcode -> Seurat cell name
  barcodes <- unique(contigs$barcode)
  mapping <- data.frame(
    barcode = barcodes,
    cell = NA_character_,
    strategy = NA_character_,
    stringsAsFactors = FALSE
  )

  # Helper to fill mapping entries if unique matches found
  fill_map <- function(idx_unset, candidates, strat) {
    if (length(idx_unset) == 0 || length(candidates) == 0) return(NULL)
    for (i in idx_unset) {
      bc <- mapping$barcode[i]
      hits <- candidates[[bc]]
      if (is.null(hits)) next
      if (length(hits) == 1L) {
        mapping$cell[i] <<- hits
        mapping$strategy[i] <<- strat
      }
    }
    invisible(NULL)
  }

  # Strategy 1: explicit cell_map
  if (!is.null(cell_map)) {
    if (is.vector(cell_map) && is.character(cell_map) && !is.null(names(cell_map))) {
      # named vector barcode -> cell
      cm <- data.frame(barcode = names(cell_map), cell = unname(cell_map), stringsAsFactors = FALSE)
    } else if (is.data.frame(cell_map) && all(c("barcode", "cell") %in% names(cell_map))) {
      cm <- data.frame(barcode = as.character(cell_map$barcode), cell = as.character(cell_map$cell), stringsAsFactors = FALSE)
    } else {
      stop("cell_map must be a named character vector (barcode->cell) or data.frame with columns 'barcode' and 'cell'.")
    }
    # Restrict to candidate cells
    cm <- cm[cm$cell %in% candidate_cells, , drop = FALSE]
    # Apply mappings where barcode is in our set
    idx <- match(cm$barcode, mapping$barcode)
    ok <- which(!is.na(idx) & !is.na(cm$cell))
    if (length(ok) > 0) {
      mapping$cell[idx[ok]] <- cm$cell[ok]
      mapping$strategy[idx[ok]] <- "cell_map"
    }
  }

  # Strategy 2: regex capture from cell names
  if (!is.null(cell_regex) && nzchar(cell_regex)) {
    extracted <- rep(NA_character_, length(candidate_cells))
    m <- regexec(cell_regex, candidate_cells)
    regmatches_list <- regmatches(candidate_cells, m)
    for (i in seq_along(regmatches_list)) {
      mm <- regmatches_list[[i]]
      if (length(mm) >= 2L) extracted[i] <- mm[2L]  # first capture group
    }
    df <- data.frame(cell = candidate_cells, barcode = extracted, stringsAsFactors = FALSE)
    df <- df[!is.na(df$barcode) & df$barcode %in% barcodes, , drop = FALSE]
    # Build candidates list: barcode -> unique cell matches
    cand <- split(df$cell, df$barcode)
    unset <- which(is.na(mapping$cell))
    fill_map(unset, cand, "cell_regex")
  }

  # Strategy 3: exact equality
  if (any(is.na(mapping$cell))) {
    hits <- intersect(barcodes, candidate_cells)
    if (length(hits) > 0) {
      idx_b <- match(hits, mapping$barcode)
      mapping$cell[idx_b] <- hits
      mapping$strategy[idx_b] <- ifelse(is.na(mapping$strategy[idx_b]), "exact", mapping$strategy[idx_b])
    }
  }

  # Strategy 4: barcode_prefix exact (prefix+barcode and prefix_+barcode)
  if (any(is.na(mapping$cell)) && !is.null(barcode_prefix) && nzchar(barcode_prefix)) {
    try1 <- paste0(barcodes, barcode_prefix) # unlikely but try both orders
    try2 <- paste0(barcode_prefix, barcodes)
    try3 <- paste0(barcode_prefix, "_", barcodes)
    cand <- lapply(barcodes, function(bc) {
      c1 <- candidate_cells[candidate_cells %in% paste0(bc, barcode_prefix)]
      c2 <- candidate_cells[candidate_cells %in% paste0(barcode_prefix, bc)]
      c3 <- candidate_cells[candidate_cells %in% paste0(barcode_prefix, "_", bc)]
      unique(c(c1, c2, c3))
    })
    names(cand) <- barcodes
    unset <- which(is.na(mapping$cell))
    fill_map(unset, cand, "barcode_prefix")
  }

  # Strategy 5: suffix match (endsWith)
  if (any(is.na(mapping$cell))) {
    # Build reverse index: for each barcode, candidate cells that end with it
    cand <- lapply(barcodes, function(bc) candidate_cells[endsWith(candidate_cells, bc)])
    names(cand) <- barcodes
    unset <- which(is.na(mapping$cell))
    fill_map(unset, cand, "suffix")
  }

  # Diagnostics based on Seurat cells queried for this run
  all_cells <- colnames(seurat)
  if (is.null(all_cells) || length(all_cells) == 0) all_cells <- rownames(seurat@meta.data)
  if (!is.null(barcode_prefix) && nzchar(barcode_prefix)) {
    queried_cells <- all_cells[startsWith(all_cells, barcode_prefix)]
  } else {
    queried_cells <- all_cells
  }
  mapped_cells <- unique(mapping$cell[!is.na(mapping$cell)])
  mapped_in_query <- intersect(queried_cells, mapped_cells)
  n_q <- length(queried_cells)
  n_matched <- length(mapped_in_query)
  n_unmatched <- n_q - n_matched
  frac_unmatched <- if (n_q > 0) n_unmatched / n_q else NA_real_
  if (!is.na(frac_unmatched) && frac_unmatched > (1 - expected_t_fraction)) {
    warning(sprintf("%d/%d queried Seurat cells did not map to VDJ clonotypes (%.1f%% unmatched).",
                    n_unmatched, n_q, 100 * frac_unmatched))
  }

  # Attach mapped cell names to contigs
  contigs$cell <- mapping$cell[match(contigs$barcode, mapping$barcode)]
  contigs$cell_id <- contigs$cell  # alias used by later steps

  # Derive per-cell, per-chain summaries and write to metadata
  md <- seurat@meta.data
  chains <- c("TRA", "TRB")
  # Helper: select best contig per (cell, chain)
  select_best_row <- function(df) {
    if (nrow(df) == 1) return(df[1, , drop = FALSE])
    # Prefer productive if present
    if ("productive" %in% names(df)) {
      prod_rows <- which(df$productive %in% c(TRUE, "True", "TRUE", 1))
      if (length(prod_rows) > 0) df <- df[prod_rows, , drop = FALSE]
    }
    # Order by UMIs desc, then Reads desc, then high_confidence if present
    umis <- if ("umis" %in% names(df)) suppressWarnings(as.numeric(df$umis)) else rep(NA_real_, nrow(df))
    reads <- if ("read_count" %in% names(df)) suppressWarnings(as.numeric(df$read_count)) else if ("reads" %in% names(df)) suppressWarnings(as.numeric(df$reads)) else rep(NA_real_, nrow(df))
    hi <- if ("high_confidence" %in% names(df)) as.integer(df$high_confidence %in% c(TRUE, "True", "TRUE", 1)) else rep(0L, nrow(df))
    ord <- order(-ifelse(is.na(umis), -Inf, umis), -ifelse(is.na(reads), -Inf, reads), -hi)
    df[ord[1], , drop = FALSE]
  }
  contigs_has <- contigs[!is.na(contigs$cell) & contigs$chain %in% chains, , drop = FALSE]
  if (nrow(contigs_has) > 0) {
    # Ensure metadata columns exist
    ensure_md_col <- function(nm, default) {
      if (!nm %in% colnames(md)) md[[nm]] <<- default
    }
    for (ch in chains) {
      ch_lower <- tolower(ch)
      sub <- contigs_has[contigs_has$chain == ch, , drop = FALSE]
      if (nrow(sub) == 0) next
      split_by_cell <- split(sub, sub$cell)
      best_list <- lapply(split_by_cell, select_best_row)
      best <- do.call(rbind, best_list)
      # Prepare values
      v_gene <- if ("v_gene" %in% names(best)) best$v_gene else NA_character_
      j_gene <- if ("j_gene" %in% names(best)) best$j_gene else NA_character_
      c_gene <- if ("c_gene" %in% names(best)) best$c_gene else NA_character_
      cdr3_aa <- if ("cdr3" %in% names(best)) best$cdr3 else if ("cdr3_aa" %in% names(best)) best$cdr3_aa else NA_character_
      cdr3_nt <- if (include_nt_subseq && "cdr3_nt" %in% names(best)) best$cdr3_nt else NULL
      umis <- if ("umis" %in% names(best)) suppressWarnings(as.numeric(best$umis)) else NA_real_
      reads <- if ("read_count" %in% names(best)) suppressWarnings(as.numeric(best$read_count)) else if ("reads" %in% names(best)) suppressWarnings(as.numeric(best$reads)) else NA_real_

      # Ensure columns exist
      ensure_md_col(paste0("vdj.", ch_lower, ".v_gene"), NA_character_)
      ensure_md_col(paste0("vdj.", ch_lower, ".j_gene"), NA_character_)
      ensure_md_col(paste0("vdj.", ch_lower, ".c_gene"), NA_character_)
      ensure_md_col(paste0("vdj.", ch_lower, ".cdr3_aa"), NA_character_)
      if (include_nt_subseq && !is.null(cdr3_nt)) ensure_md_col(paste0("vdj.", ch_lower, ".cdr3_nt"), NA_character_)
      ensure_md_col(paste0("vdj.", ch_lower, ".umis"), NA_real_)
      ensure_md_col(paste0("vdj.", ch_lower, ".reads"), NA_real_)

      rn <- rownames(best)
      md[rn, paste0("vdj.", ch_lower, ".v_gene")] <- as.character(v_gene)
      md[rn, paste0("vdj.", ch_lower, ".j_gene")] <- as.character(j_gene)
      md[rn, paste0("vdj.", ch_lower, ".c_gene")] <- as.character(c_gene)
      md[rn, paste0("vdj.", ch_lower, ".cdr3_aa")] <- as.character(cdr3_aa)
      if (include_nt_subseq && !is.null(cdr3_nt)) md[rn, paste0("vdj.", ch_lower, ".cdr3_nt")] <- as.character(cdr3_nt)
      md[rn, paste0("vdj.", ch_lower, ".umis")] <- as.numeric(umis)
      md[rn, paste0("vdj.", ch_lower, ".reads")] <- as.numeric(reads)

      # Optional CDR1/CDR2 AA and NT
      if (isTRUE(include_cdr12)) {
        if ("cdr1" %in% names(best)) {
          ensure_md_col(paste0("vdj.", ch_lower, ".cdr1_aa"), NA_character_)
          md[rn, paste0("vdj.", ch_lower, ".cdr1_aa")] <- as.character(best$cdr1)
          if (isTRUE(include_nt_subseq) && "cdr1_nt" %in% names(best)) {
            ensure_md_col(paste0("vdj.", ch_lower, ".cdr1_nt"), NA_character_)
            md[rn, paste0("vdj.", ch_lower, ".cdr1_nt")] <- as.character(best$cdr1_nt)
          }
        }
        if ("cdr2" %in% names(best)) {
          ensure_md_col(paste0("vdj.", ch_lower, ".cdr2_aa"), NA_character_)
          md[rn, paste0("vdj.", ch_lower, ".cdr2_aa")] <- as.character(best$cdr2)
          if (isTRUE(include_nt_subseq) && "cdr2_nt" %in% names(best)) {
            ensure_md_col(paste0("vdj.", ch_lower, ".cdr2_nt"), NA_character_)
            md[rn, paste0("vdj.", ch_lower, ".cdr2_nt")] <- as.character(best$cdr2_nt)
          }
        }
      }
    }
    seurat@meta.data <- md
  }

  # Prepare clonotypes with sample scoping and optional parsed sizes
  clonotypes$sample_id <- sample_id
  clonotypes$clonotype_uid <- paste0(sample_id, ":", clonotypes$clonotype_id)

  # Initialize tools[[tool_key]] if needed
  tools_list <- seurat@tools
  if (is.null(tools_list)) tools_list <- list()
  if (is.null(tools_list[[tool_key]])) tools_list[[tool_key]] <- list()

  # Append or create tables
  if (!is.null(tools_list[[tool_key]]$contigs)) {
    # Keep column union for safety
    tools_list[[tool_key]]$contigs <- rbind_fill(tools_list[[tool_key]]$contigs, contigs)
  } else {
    tools_list[[tool_key]]$contigs <- contigs
  }
  if (!is.null(tools_list[[tool_key]]$clonotypes_raw)) {
    tools_list[[tool_key]]$clonotypes_raw <- rbind_fill(tools_list[[tool_key]]$clonotypes_raw, clonotypes)
  } else {
    tools_list[[tool_key]]$clonotypes_raw <- clonotypes
  }

  # Store mapping diagnostics and configuration
  if (is.null(tools_list[[tool_key]]$mapping)) tools_list[[tool_key]]$mapping <- list()
  tools_list[[tool_key]]$mapping[[sample_id]] <- mapping
  if (is.null(tools_list[[tool_key]]$config)) tools_list[[tool_key]]$config <- list()
  tools_list[[tool_key]]$config$last_attach <- list(
    vdj_dir = vdj_dir,
    sample_id = sample_id,
    barcode_prefix = barcode_prefix,
    cell_regex = cell_regex,
    add_meta = add_meta,
    expected_t_fraction = expected_t_fraction,
    timestamp = as.character(Sys.time())
  )

  # Store per-run mapping summary
  if (is.null(tools_list[[tool_key]]$mapping_summary)) tools_list[[tool_key]]$mapping_summary <- list()
  tools_list[[tool_key]]$mapping_summary[[sample_id]] <- data.frame(
    sample_id = sample_id,
    prefix = ifelse(is.null(barcode_prefix) || !nzchar(barcode_prefix), NA_character_, barcode_prefix),
    queried = n_q,
    matched = n_matched,
    unmatched = n_unmatched,
    unmatched_fraction = frac_unmatched,
    stringsAsFactors = FALSE
  )

  # Optionally write minimal metadata: 10x clonotype id and size
  if (isTRUE(add_meta)) {
    md <- seurat@meta.data
    # Build per-barcode clonotype_id via majority across contigs; tie by summed UMIs
    has_clone <- !is.na(contigs$raw_clonotype_id) & nzchar(contigs$raw_clonotype_id)
    contigs_has <- contigs[has_clone & !is.na(contigs$cell), , drop = FALSE]
    if (nrow(contigs_has) > 0) {
      # Aggregate by cell -> clonotype_id counts and umi sums
      key <- paste(contigs_has$cell, contigs_has$raw_clonotype_id, sep = "\t")
      count_tbl <- tapply(rep(1L, length(key)), key, sum)
      umi_tbl <- tapply(as.integer(contigs_has$umis), key, sum)
      # Decide best clonotype per cell
      split_keys <- strsplit(names(count_tbl), "\t", fixed = TRUE)
      df <- data.frame(
        cell = vapply(split_keys, `[`, character(1), 1L),
        clone = vapply(split_keys, `[`, character(1), 2L),
        n = as.integer(count_tbl),
        umi = as.integer(umi_tbl[match(names(count_tbl), names(umi_tbl))]),
        stringsAsFactors = FALSE
      )
      # Order by cell, then decreasing n, then decreasing umi
      df <- df[order(df$cell, -df$n, -df$umi, df$clone), , drop = FALSE]
      best <- df[!duplicated(df$cell), c("cell", "clone")]
      # Map to sample-scoped UID and sizes
      best$clone_uid <- paste0(sample_id, ":", best$clone)
      # Lookup frequency from clonotypes table
      freq <- clonotypes$frequency
      names(freq) <- clonotypes$clonotype_uid
      best$clone_size <- as.integer(freq[best$clone_uid])

      # Ensure metadata rows exist and assign
      if (!"vdj.clone_id_10x" %in% colnames(md)) md$vdj.clone_id_10x <- NA_character_
      if (!"vdj.clone_size_10x" %in% colnames(md)) md$vdj.clone_size_10x <- NA_integer_
      rownames(best) <- best$cell
      common_cells <- intersect(rownames(md), best$cell)
      md[common_cells, "vdj.clone_id_10x"] <- best[common_cells, "clone_uid"]
      md[common_cells, "vdj.clone_size_10x"] <- best[common_cells, "clone_size"]
      seurat@meta.data <- md
    }
  }

  # Write tools slot back and return object
  seurat@tools <- tools_list
  return(seurat)
}


#' Collapse per-cell TRA/TRB into 1 TRB + up to 2 TRA pairs
#'
#' Selects one TRB per cell and up to two TRA chains ordered by support, then
#' emits a `pairs` table in `obj@tools[[tool_key]]$pairs` with one or two rows
#' per cell: the dominant pair (TRA1+TRB) and, if present and requested, the
#' secondary pair (TRA2+TRB).
#'
#' Selection heuristic:
#' - Filter contigs by QC: `productive`, `full_length`, `high_confidence` (configurable).
#' - TRB: select exactly one best contig by UMIs, then reads, then high_confidence.
#'   If more than one TRB exists and `drop_multi_trb = TRUE`, the cell is dropped.
#' - TRA: order candidate alphas by UMIs (or reads) descending, then reads, then
#'   high_confidence; keep the top one (dominant) and optionally the second.
#'
#' Metadata:
#' - Writes `vdj.tra1.*`, `vdj.tra2.*` (if present), and `vdj.trb.*` fields per cell when
#'   `write_meta = TRUE`. Also writes flags `vdj.has_tra2`, `vdj.multi_trb`, and `vdj.pair_qc`.
#' - Optionally updates the primary `vdj.tra.*` and `vdj.trb.*` fields to TRA1/TRB.
#'
#' @param seurat Seurat v5 object previously processed by `vdj_attach_10x_vdj_v2[_batch]`.
#' @param tool_key Character; tools slot name (default "vdj").
#' @param drop_multi_trb Logical; drop cells with >1 TRB after QC (default TRUE). If FALSE, picks the best TRB.
#' @param keep_two_alpha Logical; if TRUE, produce a secondary pair using TRA2+TRB when available (default TRUE).
#' @param alpha_support Character; support metric for ranking alphas: "umis" (default) or "reads".
#' @param min_umis Integer; discard contigs below this UMI count (default 0).
#' @param min_reads Integer; discard contigs below this read count (default 0).
#' @param require_productive Logical; require productive contigs (default TRUE).
#' @param require_full_length Logical; require full_length when column exists (default FALSE).
#' @param require_high_confidence Logical; require high_confidence when column exists (default FALSE).
#' @param write_meta Logical; write metadata columns described above (default TRUE).
#' @param update_primary_meta Logical; if TRUE, also set `vdj.tra.*` and `vdj.trb.*` to the primary pair (default TRUE).
#' @return Updated Seurat object (invisible).
#' @export
vdj_collapse_pairs_seurat <- function(
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
) {
  stopifnot(!is.null(seurat))
  if (!isS4(seurat) || !all(c("meta.data", "tools") %in% slotNames(seurat))) {
    stop("Expected a Seurat v5 object with slots 'meta.data' and 'tools'.")
  }
  alpha_support <- match.arg(alpha_support)

  tools_list <- seurat@tools
  if (is.null(tools_list) || is.null(tools_list[[tool_key]]) || is.null(tools_list[[tool_key]]$contigs)) {
    stop(sprintf("No contigs table found at tools$%s$contigs. Run vdj_attach_10x_vdj_v2 first.", tool_key))
  }
  contigs <- tools_list[[tool_key]]$contigs
  if (!all(c("cell", "chain", "v_gene", "j_gene", "cdr3") %in% names(contigs))) {
    stop("contigs table missing required columns: cell, chain, v_gene, j_gene, cdr3")
  }

  # Normalize QC booleans
  norm_bool <- function(x) {
    if (is.null(x) || !length(x)) return(rep(NA, nrow(contigs)))
    if (is.logical(x)) return(x)
    if (is.numeric(x)) return(x != 0)
    if (is.character(x)) return(toupper(x) %in% c("TRUE", "T", "YES", "Y"))
    as.logical(x)
  }
  contigs$productive <- if ("productive" %in% names(contigs)) norm_bool(contigs$productive) else NA
  contigs$full_length <- if ("full_length" %in% names(contigs)) norm_bool(contigs$full_length) else NA
  contigs$high_confidence <- if ("high_confidence" %in% names(contigs)) norm_bool(contigs$high_confidence) else NA
  contigs$umis <- if ("umis" %in% names(contigs)) suppressWarnings(as.numeric(contigs$umis)) else NA_real_
  contigs$reads <- if ("read_count" %in% names(contigs)) suppressWarnings(as.numeric(contigs$read_count)) else if ("reads" %in% names(contigs)) suppressWarnings(as.numeric(contigs$reads)) else NA_real_

  # Apply QC filters
  keep <- rep(TRUE, nrow(contigs))
  if (isTRUE(require_productive) && "productive" %in% names(contigs)) keep <- keep & (contigs$productive %in% c(TRUE))
  if (isTRUE(require_full_length) && "full_length" %in% names(contigs)) keep <- keep & (contigs$full_length %in% c(TRUE))
  if (isTRUE(require_high_confidence) && "high_confidence" %in% names(contigs)) keep <- keep & (contigs$high_confidence %in% c(TRUE))
  if (!is.null(min_umis) && min_umis > 0) keep <- keep & (!is.na(contigs$umis) & contigs$umis >= min_umis)
  if (!is.null(min_reads) && min_reads > 0) keep <- keep & (!is.na(contigs$reads) & contigs$reads >= min_reads)
  contigs_qc <- contigs[keep & !is.na(contigs$cell) & contigs$chain %in% c("TRA", "TRB"), , drop = FALSE]

  # Ranking helper
  rank_ord <- function(df) {
    # Prefer productive/full_length/high_confidence implicitly by filtering; in ranking use:
    umis <- if ("umis" %in% names(df)) suppressWarnings(as.numeric(df$umis)) else rep(NA_real_, nrow(df))
    reads <- if ("read_count" %in% names(df)) suppressWarnings(as.numeric(df$read_count)) else if ("reads" %in% names(df)) suppressWarnings(as.numeric(df$reads)) else rep(NA_real_, nrow(df))
    hi <- if ("high_confidence" %in% names(df)) as.integer(df$high_confidence %in% c(TRUE)) else rep(0L, nrow(df))
    if (identical(alpha_support, "umis")) {
      order(-ifelse(is.na(umis), -Inf, umis), -ifelse(is.na(reads), -Inf, reads), -hi)
    } else {
      order(-ifelse(is.na(reads), -Inf, reads), -ifelse(is.na(umis), -Inf, umis), -hi)
    }
  }

  split_by_cell <- split(contigs_qc, contigs_qc$cell)

  # Iterate cells to build pairs
  pair_rows <- list()
  md <- seurat@meta.data
  # Ensure metadata flags exist if writing
  if (isTRUE(write_meta)) {
    if (!"vdj.has_tra2" %in% colnames(md)) md$vdj.has_tra2 <- FALSE
    if (!"vdj.multi_trb" %in% colnames(md)) md$vdj.multi_trb <- FALSE
    if (!"vdj.pair_qc" %in% colnames(md)) md$vdj.pair_qc <- NA_character_
  }

  for (cell in names(split_by_cell)) {
    df <- split_by_cell[[cell]]
    tra <- df[df$chain == "TRA", , drop = FALSE]
    trb <- df[df$chain == "TRB", , drop = FALSE]

    if (nrow(trb) == 0 || nrow(tra) == 0) {
      if (isTRUE(write_meta) && cell %in% rownames(md)) md[cell, "vdj.pair_qc"] <- if (nrow(trb) == 0) "no_trb" else "no_tra"
      next
    }

    # Handle multiple TRB
    multi_trb <- nrow(trb) > 1
    if (isTRUE(write_meta) && cell %in% rownames(md)) md[cell, "vdj.multi_trb"] <- multi_trb
    if (multi_trb && isTRUE(drop_multi_trb)) {
      if (isTRUE(write_meta) && cell %in% rownames(md)) md[cell, "vdj.pair_qc"] <- "multi_trb_dropped"
      next
    }

    # Pick best TRB (even if multiple and not dropped)
    trb <- trb[rank_ord(trb), , drop = FALSE]
    trb_best <- trb[1, , drop = FALSE]

    # Order TRA by support and pick 1 or 2
    tra <- tra[rank_ord(tra), , drop = FALSE]
    tra_best <- tra[1, , drop = FALSE]
    tra_second <- if (isTRUE(keep_two_alpha) && nrow(tra) >= 2) tra[2, , drop = FALSE] else NULL

    # Build primary pair row
    mk_row <- function(rank, tra_row, trb_row) {
      tra_reads_val <- NA
      trb_reads_val <- NA
      if ("read_count" %in% names(tra_row)) tra_reads_val <- suppressWarnings(as.numeric(tra_row$read_count))
      else if ("reads" %in% names(tra_row)) tra_reads_val <- suppressWarnings(as.numeric(tra_row$reads))
      if ("read_count" %in% names(trb_row)) trb_reads_val <- suppressWarnings(as.numeric(trb_row$read_count))
      else if ("reads" %in% names(trb_row)) trb_reads_val <- suppressWarnings(as.numeric(trb_row$reads))
      data.frame(
        cell = cell,
        pair_rank = as.integer(rank),
        tra_cdr3 = as.character(tra_row$cdr3),
        tra_v = as.character(tra_row$v_gene),
        tra_j = as.character(tra_row$j_gene),
        tra_umis = suppressWarnings(as.numeric(tra_row$umis)),
        tra_reads = tra_reads_val,
        trb_cdr3 = as.character(trb_row$cdr3),
        trb_v = as.character(trb_row$v_gene),
        trb_j = as.character(trb_row$j_gene),
        trb_umis = suppressWarnings(as.numeric(trb_row$umis)),
        trb_reads = trb_reads_val,
        stringsAsFactors = FALSE
      )
    }
    pair_rows[[length(pair_rows) + 1L]] <- mk_row(1L, tra_best, trb_best)
    if (!is.null(tra_second)) pair_rows[[length(pair_rows) + 1L]] <- mk_row(2L, tra_second, trb_best)

    if (isTRUE(write_meta) && cell %in% rownames(md)) {
      # Ensure metadata columns
      ensure_md <- function(nm, default) { if (!nm %in% colnames(md)) md[[nm]] <<- default }
      # Primary
      ensure_md("vdj.tra1.cdr3_aa", NA_character_)
      ensure_md("vdj.tra1.v_gene", NA_character_)
      ensure_md("vdj.tra1.j_gene", NA_character_)
      ensure_md("vdj.trb.cdr3_aa", NA_character_)
      ensure_md("vdj.trb.v_gene", NA_character_)
      ensure_md("vdj.trb.j_gene", NA_character_)
      md[cell, "vdj.tra1.cdr3_aa"] <- as.character(tra_best$cdr3)
      md[cell, "vdj.tra1.v_gene"] <- as.character(tra_best$v_gene)
      md[cell, "vdj.tra1.j_gene"] <- as.character(tra_best$j_gene)
      md[cell, "vdj.trb.cdr3_aa"] <- as.character(trb_best$cdr3)
      md[cell, "vdj.trb.v_gene"] <- as.character(trb_best$v_gene)
      md[cell, "vdj.trb.j_gene"] <- as.character(trb_best$j_gene)
      # Secondary
      ensure_md("vdj.tra2.cdr3_aa", NA_character_)
      ensure_md("vdj.tra2.v_gene", NA_character_)
      ensure_md("vdj.tra2.j_gene", NA_character_)
      md[cell, "vdj.has_tra2"] <- !is.null(tra_second)
      if (!is.null(tra_second)) {
        md[cell, "vdj.tra2.cdr3_aa"] <- as.character(tra_second$cdr3)
        md[cell, "vdj.tra2.v_gene"] <- as.character(tra_second$v_gene)
        md[cell, "vdj.tra2.j_gene"] <- as.character(tra_second$j_gene)
      } else {
        md[cell, c("vdj.tra2.cdr3_aa", "vdj.tra2.v_gene", "vdj.tra2.j_gene")] <- NA_character_
      }
      md[cell, "vdj.pair_qc"] <- "pass"

      if (isTRUE(update_primary_meta)) {
        # Also maintain legacy fields to primary pair
        ensure_md("vdj.tra.cdr3_aa", NA_character_)
        ensure_md("vdj.tra.v_gene", NA_character_)
        ensure_md("vdj.tra.j_gene", NA_character_)
        md[cell, "vdj.tra.cdr3_aa"] <- md[cell, "vdj.tra1.cdr3_aa"]
        md[cell, "vdj.tra.v_gene"] <- md[cell, "vdj.tra1.v_gene"]
        md[cell, "vdj.tra.j_gene"] <- md[cell, "vdj.tra1.j_gene"]
      }
    }
  }

  pairs <- if (length(pair_rows)) do.call(rbind, pair_rows) else data.frame()

  # Store pairs table
  tools_list[[tool_key]]$pairs <- pairs
  seurat@tools <- tools_list

  # Write metadata if changed
  if (isTRUE(write_meta)) seurat@meta.data <- md
  invisible(seurat)
}

# Internal: rbind with column union using base R
rbind_fill <- function(x, y) {
  if (is.null(x) || nrow(x) == 0) return(y)
  if (is.null(y) || nrow(y) == 0) return(x)
  xn <- names(x); yn <- names(y)
  add_to_x <- setdiff(yn, xn)
  add_to_y <- setdiff(xn, yn)
  if (length(add_to_x) > 0) for (nm in add_to_x) x[[nm]] <- NA
  if (length(add_to_y) > 0) for (nm in add_to_y) y[[nm]] <- NA
  # Reorder columns identically
  y <- y[, names(x), drop = FALSE]
  rbind(x, y)
}


#' Attach multiple 10x 5' VDJ v2 runs using a caps-token template and prefix map
#'
#' Convenience wrapper over `vdj_attach_10x_vdj_v2()` to batch-attach many
#' Cell Ranger VDJ outputs to a single Seurat object when cell names were
#' prefixed to avoid barcode collisions. The user supplies a data.frame with a
#' prefix column and any number of UPPERCASE columns that correspond to tokens in
#' a path template. All UPPERCASE tokens present in the template (e.g.,
#' `FOLDERNAME1`, `FOLDERNAME2`) are replaced per-row by the values of columns of
#' the same names to locate each run's VDJ output directory.
#'
#' Example template with two tokens:
#' `/path/to/MMCarT_240212_reseq/FOLDERNAME1/data/FOLDERNAME2/outs/per_sample_outs/FOLDERNAME2/vdj_t`
#'
#' If no UPPERCASE tokens are detected in the template, the function falls back
#' to single-placeholder mode, replacing `placeholder` with values from
#' `folders_col`.
#'
#' @param seurat A Seurat v5 object.
#' @param prefix_folders data.frame with at least a prefix column and one column per UPPERCASE token present in `vdj_dir_template`.
#' @param vdj_dir_template Character scalar path containing UPPERCASE token(s) such as `FOLDERNAME1`, `FOLDERNAME2` to be substituted from `prefix_folders` per row.
#' @param tool_key Name under `obj@tools` to store results (default `"vdj"`).
#' @param add_meta If TRUE, each run writes `vdj.clone_id_10x` and `vdj.clone_size_10x` for mapped cells.
#' @param include_cdr12 If TRUE, forward to per-run attach to import CDR1/2 AA where available.
#' @param include_nt_subseq If TRUE, forward to per-run attach to include nucleotide sequences for CDR1/2/3 where available.
#' @param expected_t_fraction Expected fraction of T cells among queried Seurat cells; forwarded to per-run attach.
#'   A warning is emitted only if `(unmatched / queried) > (1 - expected_t_fraction)`. Default: `0.5`.
#' @param placeholder Token to replace in the template when no UPPERCASE tokens are detected (default `"FOLDERNAME"`).
#' @param prefixes_col Column name(s) in `prefix_folders` to look for the prefix. First found is used (default tries `"prefixes"`, then `"PREFIXES"`).
#' @param folders_col Column name in `prefix_folders` for folder names used only in single-placeholder fallback mode (default `"folders"`).
#' @param cell_regex Optional regex with one capture group to extract barcodes from cell names (forwarded to per-run attach).
#' @param cell_map Optional explicit mapping (forwarded to per-run attach).
#'
#' @return Updated Seurat object (invisible).
#' @export
vdj_attach_10x_vdj_v2_batch <- function(
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
) {
  stopifnot(!is.null(seurat))
  if (!isS4(seurat) || !all(c("meta.data", "tools") %in% slotNames(seurat))) {
    stop("Expected a Seurat v5 object with slots 'meta.data' and 'tools'.")
  }

  # Validate mapping frame and template
  if (!is.data.frame(prefix_folders)) {
    stop("prefix_folders must be a data.frame with prefix and folder columns.")
  }
  # Resolve prefixes column with fallback
  pf_cols <- colnames(prefix_folders)
  pref_col <- NULL
  if (length(prefixes_col) > 1L) {
    for (cand in prefixes_col) {
      if (cand %in% pf_cols) { pref_col <- cand; break }
    }
  } else if (length(prefixes_col) == 1L && prefixes_col %in% pf_cols) {
    pref_col <- prefixes_col
  }
  if (is.null(pref_col)) {
    stop(sprintf(
      "prefix_folders must contain a prefix column; looked for: %s",
      paste(prefixes_col, collapse = ", ")
    ))
  }
  if (!is.character(vdj_dir_template) || length(vdj_dir_template) != 1L || !nzchar(vdj_dir_template)) {
    stop("vdj_dir_template must be a non-empty character scalar containing the placeholder token.")
  }
  # Detect CAPPED tokens based on columns actually present in the mapping
  caps_cols <- pf_cols[grepl("^[A-Z0-9_]+$", pf_cols)]
  caps_tokens <- caps_cols[vapply(caps_cols, function(nm) grepl(nm, vdj_dir_template, fixed = TRUE), logical(1))]
  # Prefer longer tokens first to avoid partial overlaps
  caps_tokens <- caps_tokens[order(-nchar(caps_tokens))]
  # If no caps tokens found, fall back to single placeholder behavior
  use_caps <- length(caps_tokens) > 0L
  if (!use_caps) {
    if (!nzchar(placeholder) || !grepl(placeholder, vdj_dir_template, fixed = TRUE)) {
      stop(sprintf("Template must contain placeholder token '%s' at least once.", placeholder))
    }
    if (!folders_col %in% pf_cols) {
      stop(sprintf("prefix_folders is missing required column '%s' for placeholder substitution.", folders_col))
    }
  }

  # Iterate rows and attach
  n <- nrow(prefix_folders)
  if (n == 0) return(seurat)
  for (i in seq_len(n)) {
    prefix <- as.character(prefix_folders[[pref_col]][i])
    if (!nzchar(prefix)) next
    if (use_caps) {
      vdj_dir <- vdj_dir_template
      for (tok in caps_tokens) {
        val <- as.character(prefix_folders[[tok]][i])
        if (!nzchar(val)) {
          warning(sprintf("Skipping row %d: empty value for token '%s'", i, tok))
          vdj_dir <- NA_character_
          break
        }
        vdj_dir <- gsub(tok, val, vdj_dir, fixed = TRUE)
      }
      if (is.na(vdj_dir)) next
    } else {
      folder <- as.character(prefix_folders[[folders_col]][i])
      if (!nzchar(folder)) next
      vdj_dir <- gsub(placeholder, folder, vdj_dir_template, fixed = TRUE)
    }
    # Skip if directory missing but warn for visibility
    if (!dir.exists(vdj_dir)) {
      if (use_caps) {
        vals <- vapply(caps_tokens, function(tk) as.character(prefix_folders[[tk]][i]), character(1))
        msg_kv <- paste(sprintf("%s=%s", caps_tokens, vals), collapse = ", ")
        warning(sprintf("Skipping row %d: VDJ directory not found (%s) -> %s", i, msg_kv, vdj_dir))
      } else {
        folder <- as.character(prefix_folders[[folders_col]][i])
        warning(sprintf("Skipping: VDJ directory not found for folder '%s' -> %s", folder, vdj_dir))
      }
      next
    }
    # Delegate to single-run attach; let it infer sample_id from the path
    seurat <- vdj_attach_10x_vdj_v2(
      seurat = seurat,
      vdj_dir = vdj_dir,
      sample_id = NULL,
      barcode_prefix = prefix,
      cell_regex = cell_regex,
      cell_map = cell_map,
      tool_key = tool_key,
      add_meta = add_meta,
      include_cdr12 = include_cdr12,
      include_nt_subseq = include_nt_subseq,
      expected_t_fraction = expected_t_fraction
    )
  }
  return(seurat)
}
