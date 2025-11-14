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
#' @param progress show progress bar (default TRUE)
#' @param chunk_size number of queries to process per chunk (default 5000)
#' @return data.frame with query metadata and hit columns
#' @export
match_tcr_many_df <- function(db, cdr3, v_segment, j_segment, scope = "0,0,0,0", top_n = 0L,
                               progress = TRUE, chunk_size = 5000L) {
  n_queries <- length(cdr3)

  # For small batches, just run directly without chunking
  if (n_queries <= chunk_size || !progress) {
    res <- match_tcr_many(db, as.character(cdr3), as.character(v_segment),
                          as.character(j_segment), scope, as.integer(top_n))
    return(as.data.frame(res, stringsAsFactors = FALSE))
  }

  # Chunk processing with progress bar
  n_chunks <- ceiling(n_queries / chunk_size)

  if (progress) {
    message(sprintf("Matching %d TCRs against database (parallel processing enabled)", n_queries))
    pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  }

  results_list <- list()

  for (i in seq_len(n_chunks)) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_queries)
    idx <- start_idx:end_idx

    # Match this chunk (uses parallel processing in Rust)
    chunk_res <- match_tcr_many(
      db,
      as.character(cdr3[idx]),
      as.character(v_segment[idx]),
      as.character(j_segment[idx]),
      scope,
      as.integer(top_n)
    )

    chunk_df <- as.data.frame(chunk_res, stringsAsFactors = FALSE)

    # Adjust query_index to global indices
    if (nrow(chunk_df) > 0 && "query_index" %in% colnames(chunk_df)) {
      chunk_df$query_index <- chunk_df$query_index + start_idx - 1
    }

    results_list[[i]] <- chunk_df

    if (progress) setTxtProgressBar(pb, i)
  }

  if (progress) {
    close(pb)
    message(sprintf("Matched %d queries, found %d hits",
                    n_queries, sum(sapply(results_list, nrow))))
  }

  # Combine all chunks
  do.call(rbind, results_list)
}
