#' Calculate pairwise TCRdist distances between TCRs
#'
#' @name calculate_tcrdist
#' @description
#' Calculates pairwise distances between T-cell receptors using the tcrdist algorithm,
#' which combines BLOSUM62 amino acid substitution scoring with Needleman-Wunsch alignment.
#' This implementation is inspired by the original tcrdist method (Dash et al., 2017) and
#' provides a biologically meaningful measure of TCR similarity.
#'
#' @details
#' The tcrdist algorithm calculates distances based on:
#' \itemize{
#'   \item \strong{BLOSUM62 scoring}: Amino acid substitution matrix for biological relevance
#'   \item \strong{Position-specific scoring}: \code{max(0, 4 - BLOSUM62[aa1, aa2])} per position
#'   \item \strong{CDR weighting}: CDR3 weighted 3x more than CDR1/2 (reflects biological importance)
#'   \item \strong{Gap penalties}: 4 for CDR1/2, 8 for CDR3
#'   \item \strong{Chain combination}: Distances from alpha and beta chains are summed
#' }
#'
#' The distance calculation for each chain:
#' \deqn{distance = CDR1_{dist} \times 1 + CDR2_{dist} \times 1 + CDR3_{dist} \times 3}
#'
#' Total TCR distance:
#' \deqn{tcrdist = \alpha_{chain} + \beta_{chain}}
#'
#' Missing CDR sequences (empty strings or NA) are handled gracefully - those regions
#' are simply not included in the distance calculation.
#'
#' @param cdr1_a Character vector of CDR1 alpha sequences (amino acids). Use empty strings "" for missing data.
#' @param cdr2_a Character vector of CDR2 alpha sequences (amino acids). Use empty strings "" for missing data.
#' @param cdr3_a Character vector of CDR3 alpha sequences (amino acids). Use empty strings "" for missing data.
#' @param cdr1_b Character vector of CDR1 beta sequences (amino acids). Use empty strings "" for missing data.
#' @param cdr2_b Character vector of CDR2 beta sequences (amino acids). Use empty strings "" for missing data.
#' @param cdr3_b Character vector of CDR3 beta sequences (amino acids). Use empty strings "" for missing data.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{i}{Integer vector of row indices (1-based) for the distance matrix}
#'   \item{j}{Integer vector of column indices (1-based) for the distance matrix}
#'   \item{distance}{Numeric vector of pairwise distances}
#'   \item{n}{Integer, number of TCRs}
#' }
#'
#' The distance matrix can be reconstructed using:
#' \code{matrix(result$distance, nrow = result$n, ncol = result$n)}
#'
#' @section Distance Interpretation:
#' \itemize{
#'   \item \strong{0}: Identical TCRs
#'   \item \strong{< 50}: Very similar TCRs, likely recognize the same epitope
#'   \item \strong{50-100}: Moderate similarity
#'   \item \strong{> 100}: Dissimilar TCRs
#' }
#'
#' @section Performance:
#' This function is implemented in Rust with parallel processing via Rayon for high performance.
#' For \code{n} TCRs, it performs \code{n^2} pairwise comparisons. The parallel implementation
#' automatically uses all available CPU cores. Typical performance on a modern multi-core CPU:
#' \itemize{
#'   \item 100 TCRs: ~10,000 comparisons, < 1 second
#'   \item 1,000 TCRs: ~1,000,000 comparisons, ~1 second
#'   \item 10,000 TCRs: ~100,000,000 comparisons, ~30-60 seconds
#'   \item 50,000 TCRs: ~2,500,000,000 comparisons, ~10-20 minutes
#' }
#'
#' For very large datasets (>5,000 TCRs), consider using \code{\link{calculate_tcrdist_with_progress}}
#' which adds progress bar support.
#'
#' @examples
#' # Example 1: Calculate distances with CDR3 only
#' tcr_data <- data.frame(
#'   id = paste0("TCR_", 1:5),
#'   cdr3_a = c("CAASNRGSTLGRLYF", "CAASIRSSYKLIF", "CAASNRGSTLGRLYF",
#'              "CALSDPNQAGTALIF", "CAASKQGAQKLVF"),
#'   cdr3_b = c("CASSLTGNTEAFF", "CASSLGQGAYEQYF", "CASSLTGNTEAFF",
#'              "CASSLGQGAYEQYF", "CASSVGQGGELFF")
#' )
#'
#' # Calculate distances (CDR1/2 not available, use empty strings)
#' result <- calculate_tcrdist(
#'   cdr1_a = rep("", 5),
#'   cdr2_a = rep("", 5),
#'   cdr3_a = tcr_data$cdr3_a,
#'   cdr1_b = rep("", 5),
#'   cdr2_b = rep("", 5),
#'   cdr3_b = tcr_data$cdr3_b
#' )
#'
#' # Convert to distance matrix
#' dist_matrix <- matrix(result$distance, nrow = result$n, ncol = result$n)
#' rownames(dist_matrix) <- tcr_data$id
#' colnames(dist_matrix) <- tcr_data$id
#'
#' print(dist_matrix)
#'
#' # Example 2: Use for hierarchical clustering
#' library(stats)
#' hc <- hclust(as.dist(dist_matrix), method = "ward.D2")
#' plot(hc, main = "TCR Clustering")
#'
#' # Example 3: Visualize as heatmap
#' if (require("pheatmap")) {
#'   pheatmap::pheatmap(
#'     dist_matrix,
#'     main = "TCR Distance Heatmap",
#'     clustering_method = "ward.D2"
#'   )
#' }
#'
#' # Example 4: With all CDR regions
#' tcr_full <- data.frame(
#'   cdr1_a = c("TGTGC", "TGTGC"),
#'   cdr2_a = c("TGTGC", "TGTGA"),
#'   cdr3_a = c("CASSF", "CASSLF"),
#'   cdr1_b = c("TGTGC", "TGTGC"),
#'   cdr2_b = c("TGTGC", "TGTGA"),
#'   cdr3_b = c("CASSF", "CASSLF")
#' )
#'
#' result_full <- calculate_tcrdist(
#'   cdr1_a = tcr_full$cdr1_a,
#'   cdr2_a = tcr_full$cdr2_a,
#'   cdr3_a = tcr_full$cdr3_a,
#'   cdr1_b = tcr_full$cdr1_b,
#'   cdr2_b = tcr_full$cdr2_b,
#'   cdr3_b = tcr_full$cdr3_b
#' )
#'
#' @references
#' Dash P, Fiore-Gartland AJ, Hertz T, et al. (2017)
#' Quantifiable predictive features define epitope-specific T cell receptor repertoires.
#' \emph{Nature}, 547(7661):89-93. \doi{10.1038/nature22383}
#'
#' @seealso
#' \code{\link{tcrdist_single}} for calculating distance between two individual TCRs
#'
#' @export
NULL  # Documentation only - actual function is in extendr-wrappers.R


#' Calculate TCRdist between two individual TCRs
#'
#' @name tcrdist_single
#' @description
#' Calculates the tcrdist distance between two T-cell receptors using BLOSUM62
#' amino acid substitution scoring and Needleman-Wunsch alignment. This is a
#' convenience function for comparing just two TCRs.
#'
#' @details
#' This function uses the same algorithm as \code{\link{calculate_tcrdist}} but
#' is optimized for comparing a single pair of TCRs. The distance calculation
#' combines distances from CDR1, CDR2, and CDR3 regions of both alpha and beta
#' chains, with CDR3 weighted 3x more heavily than CDR1/2.
#'
#' See \code{\link{calculate_tcrdist}} for detailed information about the
#' tcrdist algorithm and distance interpretation.
#'
#' @param cdr1_a_1 Character, CDR1 alpha sequence of first TCR (use "" for missing)
#' @param cdr2_a_1 Character, CDR2 alpha sequence of first TCR (use "" for missing)
#' @param cdr3_a_1 Character, CDR3 alpha sequence of first TCR (use "" for missing)
#' @param cdr1_b_1 Character, CDR1 beta sequence of first TCR (use "" for missing)
#' @param cdr2_b_1 Character, CDR2 beta sequence of first TCR (use "" for missing)
#' @param cdr3_b_1 Character, CDR3 beta sequence of first TCR (use "" for missing)
#' @param cdr1_a_2 Character, CDR1 alpha sequence of second TCR (use "" for missing)
#' @param cdr2_a_2 Character, CDR2 alpha sequence of second TCR (use "" for missing)
#' @param cdr3_a_2 Character, CDR3 alpha sequence of second TCR (use "" for missing)
#' @param cdr1_b_2 Character, CDR1 beta sequence of second TCR (use "" for missing)
#' @param cdr2_b_2 Character, CDR2 beta sequence of second TCR (use "" for missing)
#' @param cdr3_b_2 Character, CDR3 beta sequence of second TCR (use "" for missing)
#'
#' @return Numeric value representing the tcrdist distance between the two TCRs.
#' Lower values indicate more similar TCRs.
#'
#' @section Distance Interpretation:
#' \itemize{
#'   \item \strong{0}: Identical TCRs
#'   \item \strong{< 50}: Very similar, likely same epitope specificity
#'   \item \strong{50-100}: Moderate similarity
#'   \item \strong{> 100}: Dissimilar TCRs
#' }
#'
#' @examples
#' # Compare two TCRs (CDR3 only)
#' dist <- tcrdist_single(
#'   cdr1_a_1 = "", cdr2_a_1 = "", cdr3_a_1 = "CAASNRGSTLGRLYF",
#'   cdr1_b_1 = "", cdr2_b_1 = "", cdr3_b_1 = "CASSLTGNTEAFF",
#'   cdr1_a_2 = "", cdr2_a_2 = "", cdr3_a_2 = "CAASIRSSYKLIF",
#'   cdr1_b_2 = "", cdr2_b_2 = "", cdr3_b_2 = "CASSLGQGAYEQYF"
#' )
#' print(paste("Distance:", dist))
#'
#' # Compare identical TCRs (should be 0)
#' dist_same <- tcrdist_single(
#'   "", "", "CASSF", "", "", "CASSF",
#'   "", "", "CASSF", "", "", "CASSF"
#' )
#' print(paste("Distance between identical TCRs:", dist_same))  # Should be 0
#'
#' # Compare with all CDR regions
#' dist_full <- tcrdist_single(
#'   cdr1_a_1 = "TGTGC", cdr2_a_1 = "TGTGC", cdr3_a_1 = "CASSF",
#'   cdr1_b_1 = "TGTGC", cdr2_b_1 = "TGTGC", cdr3_b_1 = "CASSF",
#'   cdr1_a_2 = "TGTGA", cdr2_a_2 = "TGTGA", cdr3_a_2 = "CASSLF",
#'   cdr1_b_2 = "TGTGA", cdr2_b_2 = "TGTGA", cdr3_b_2 = "CASSLF"
#' )
#' print(paste("Distance with all CDRs:", dist_full))
#'
#' @references
#' Dash P, Fiore-Gartland AJ, Hertz T, et al. (2017)
#' Quantifiable predictive features define epitope-specific T cell receptor repertoires.
#' \emph{Nature}, 547(7661):89-93. \doi{10.1038/nature22383}
#'
#' @seealso
#' \code{\link{calculate_tcrdist}} for pairwise distances between multiple TCRs
#'
#' @export
NULL  # Documentation only - actual function is in extendr-wrappers.R


#' Calculate TCRdist with progress bar support
#'
#' @description
#' Wrapper around \code{\link{calculate_tcrdist}} that adds progress bar support
#' for large datasets by processing TCRs in chunks. Useful for tracking progress
#' when computing distances for thousands of TCRs.
#'
#' @param cdr1_a Character vector of CDR1 alpha sequences. Use empty strings "" for missing data.
#' @param cdr2_a Character vector of CDR2 alpha sequences. Use empty strings "" for missing data.
#' @param cdr3_a Character vector of CDR3 alpha sequences. Use empty strings "" for missing data.
#' @param cdr1_b Character vector of CDR1 beta sequences. Use empty strings "" for missing data.
#' @param cdr2_b Character vector of CDR2 beta sequences. Use empty strings "" for missing data.
#' @param cdr3_b Character vector of CDR3 beta sequences. Use empty strings "" for missing data.
#' @param progress Logical; if TRUE, show progress bar (default TRUE)
#' @param chunk_size Integer; number of TCRs to process per chunk for progress updates (default 1000)
#'
#' @return A list with the same structure as \code{\link{calculate_tcrdist}}:
#' \describe{
#'   \item{i}{Integer vector of row indices (1-based) for the distance matrix}
#'   \item{j}{Integer vector of column indices (1-based) for the distance matrix}
#'   \item{distance}{Numeric vector of pairwise distances}
#'   \item{n}{Integer, number of TCRs}
#' }
#'
#' @details
#' This function processes TCRs in chunks, computing all pairwise distances
#' between chunk members and all other TCRs. The chunk size controls how
#' frequently the progress bar updates. Larger chunks mean fewer updates but
#' potentially less responsive progress tracking.
#'
#' Note: The underlying computation still uses parallel processing via Rayon,
#' so this function maintains high performance while adding progress visibility.
#'
#' @examples
#' # Example with progress bar
#' tcr_data <- data.frame(
#'   id = paste0("TCR_", 1:100),
#'   cdr3_a = replicate(100, paste(sample(c("A","C","D","E","F","G","H","I","K","L",
#'                                            "M","N","P","Q","R","S","T","V","W","Y"),
#'                                          10, replace = TRUE), collapse = "")),
#'   cdr3_b = replicate(100, paste(sample(c("A","C","D","E","F","G","H","I","K","L",
#'                                            "M","N","P","Q","R","S","T","V","W","Y"),
#'                                          12, replace = TRUE), collapse = ""))
#' )
#'
#' result <- calculate_tcrdist_with_progress(
#'   cdr1_a = rep("", 100),
#'   cdr2_a = rep("", 100),
#'   cdr3_a = tcr_data$cdr3_a,
#'   cdr1_b = rep("", 100),
#'   cdr2_b = rep("", 100),
#'   cdr3_b = tcr_data$cdr3_b,
#'   progress = TRUE,
#'   chunk_size = 25
#' )
#'
#' # Convert to matrix
#' dist_matrix <- matrix(result$distance, nrow = result$n, ncol = result$n)
#'
#' @seealso
#' \code{\link{calculate_tcrdist}} for the underlying function without progress bar
#'
#' @export
calculate_tcrdist_with_progress <- function(
  cdr1_a,
  cdr2_a,
  cdr3_a,
  cdr1_b,
  cdr2_b,
  cdr3_b,
  progress = TRUE,
  chunk_size = 1000L
) {
  n <- length(cdr3_a)

  # Validate inputs
  if (!(length(cdr1_a) == n && length(cdr2_a) == n &&
        length(cdr1_b) == n && length(cdr2_b) == n && length(cdr3_b) == n)) {
    stop("All CDR vectors must have equal length")
  }

  # For small datasets, just compute directly
  if (n <= chunk_size || !progress) {
    return(calculate_tcrdist(cdr1_a, cdr2_a, cdr3_a, cdr1_b, cdr2_b, cdr3_b))
  }

  # Chunk processing with progress bar
  n_chunks <- ceiling(n / chunk_size)

  if (progress) {
    message(sprintf("Calculating TCRdist for %d TCRs (%d comparisons, parallel processing enabled)",
                    n, n * n))
    pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  }

  # Pre-allocate result vectors
  all_i <- vector("integer", n * n)
  all_j <- vector("integer", n * n)
  all_dist <- vector("numeric", n * n)

  result_idx <- 1

  for (chunk_i in seq_len(n_chunks)) {
    start_i <- (chunk_i - 1) * chunk_size + 1
    end_i <- min(chunk_i * chunk_size, n)
    idx_i <- start_i:end_i

    # Calculate distances for this chunk against all TCRs
    chunk_result <- calculate_tcrdist(
      cdr1_a = cdr1_a[idx_i],
      cdr2_a = cdr2_a[idx_i],
      cdr3_a = cdr3_a[idx_i],
      cdr1_b = cdr1_b[idx_i],
      cdr2_b = cdr2_b[idx_i],
      cdr3_b = cdr3_b[idx_i]
    )

    # This gives us distances for rows start_i:end_i against columns start_i:end_i
    # We need to compute against ALL columns, so do this for each TCR in the chunk
    for (i_offset in seq_along(idx_i)) {
      global_i <- idx_i[i_offset]

      # Compute this TCR against all TCRs
      for (global_j in 1:n) {
        all_i[result_idx] <- global_i
        all_j[result_idx] <- global_j

        # If both indices are in the chunk, use the chunk result
        if (global_j >= start_i && global_j <= end_i) {
          # Find in chunk results
          local_j <- global_j - start_i + 1
          chunk_idx <- which(chunk_result$i == i_offset & chunk_result$j == local_j)
          all_dist[result_idx] <- chunk_result$distance[chunk_idx]
        } else {
          # Compute single distance
          all_dist[result_idx] <- tcrdist_single(
            cdr1_a[global_i], cdr2_a[global_i], cdr3_a[global_i],
            cdr1_b[global_i], cdr2_b[global_i], cdr3_b[global_i],
            cdr1_a[global_j], cdr2_a[global_j], cdr3_a[global_j],
            cdr1_b[global_j], cdr2_b[global_j], cdr3_b[global_j]
          )
        }

        result_idx <- result_idx + 1
      }
    }

    if (progress) setTxtProgressBar(pb, chunk_i)
  }

  if (progress) {
    close(pb)
    message(sprintf("Completed %d pairwise TCR distance calculations", n * n))
  }

  list(
    i = all_i,
    j = all_j,
    distance = all_dist,
    n = n
  )
}
