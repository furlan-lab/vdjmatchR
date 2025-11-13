#' Match one clonotype and return a data.frame
#'
#' @param db an RDatabase object
#' @param cdr3 CDR3 amino-acid sequence
#' @param v_segment V segment (optional; empty string to ignore)
#' @param j_segment J segment (optional; empty string to ignore)
#' @param scope search scope string like "0,0,0,0" or "2,1,2,3"
#' @param top_n keep top N hits (per query)
#' @return data.frame with matching hits
#' @export
match_tcr_df <- function(db, cdr3, v_segment = "", j_segment = "", scope = "0,0,0,0", top_n = 0L) {
  res <- match_tcr(db, cdr3, v_segment, j_segment, scope, as.integer(top_n))
  as.data.frame(res, stringsAsFactors = FALSE)
}

#' Match many clonotypes and return a data.frame stacked across queries
#'
#' @param db an RDatabase object
#' @param cdr3 character vector of CDR3 sequences
#' @param v_segment character vector of V segments (same length)
#' @param j_segment character vector of J segments (same length)
#' @param scope search scope string like "0,0,0,0" or "2,1,2,3"
#' @param top_n keep top N hits per query
#' @return data.frame with query metadata and hit columns
#' @export
match_tcr_many_df <- function(db, cdr3, v_segment, j_segment, scope = "0,0,0,0", top_n = 0L) {
  res <- match_tcr_many(db, as.character(cdr3), as.character(v_segment), as.character(j_segment), scope, as.integer(top_n))
  as.data.frame(res, stringsAsFactors = FALSE)
}
