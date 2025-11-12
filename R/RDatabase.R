#' RDatabase: Handle to a loaded VDJdb
#'
#' @description
#' `RDatabase` is an object provided by the Rust backend (via extendr) that
#' holds an in-memory VDJdb table and exposes methods for filtering.
#'
#' Methods include:
#' - `new_from_file(path)`
#' - `new_from_vdjdb(use_fat_db = FALSE)`
#' - `len()`
#' - `filter(species = NULL, gene = NULL, min_vdjdb_score = 0)`
#' - `filter_by_epitope_size(min_size = 1)`
#'
#' See also: [match_tcr_df()], [match_tcr_many_df()], [vdjdb_download()].
#'
#' @name RDatabase
#' @docType class
#' @keywords internal
NULL
