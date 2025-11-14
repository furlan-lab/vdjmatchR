#' Convert VDJdb database to data.table
#'
#' Extracts all entries from the database and returns as a data.table for easy
#' inspection and manipulation in R.
#'
#' @param db an RDatabase object
#' @return data.table with database entries
#' @export
#' @examples
#' \dontrun{
#' # Load database
#' db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))
#'
#' # Convert to data.table
#' dt <- db_to_table(db)
#'
#' # Filter human TRB entries
#' db_filtered <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 0)
#' dt_trb <- db_to_table(db_filtered)
#'
#' # Now use data.table operations
#' library(data.table)
#' dt_trb[, .N, by = antigen_epitope][order(-N)][1:10]  # Top 10 epitopes
#' }
db_to_table <- function(db) {
  if (!inherits(db, "RDatabase")) {
    stop("db must be an RDatabase object (created with vdjdb_open_file)")
  }

  # Get columns from Rust
  cols <- db$to_columns()

  # Convert to data.table if available, otherwise data.frame
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt <- data.table::setDT(as.data.frame(cols, stringsAsFactors = FALSE))
    return(dt)
  } else {
    # Fall back to data.frame
    df <- as.data.frame(cols, stringsAsFactors = FALSE)
    message("Install data.table for better performance: install.packages('data.table')")
    return(df)
  }
}

#' Convert VDJdb database to data.frame
#'
#' Extracts all entries from the database and returns as a data.frame.
#' For larger databases, consider using \code{db_to_table()} which returns
#' a data.table for better performance.
#'
#' @param db an RDatabase object
#' @return data.frame with database entries
#' @export
#' @examples
#' \dontrun{
#' db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = FALSE))
#' df <- db_to_df(db)
#' head(df)
#' }
db_to_df <- function(db) {
  if (!inherits(db, "RDatabase")) {
    stop("db must be an RDatabase object (created with vdjdb_open_file)")
  }

  cols <- db$to_columns()
  as.data.frame(cols, stringsAsFactors = FALSE)
}

#' Quick database summary
#'
#' Print a summary of database contents by species and gene
#'
#' @param db an RDatabase object
#' @return invisibly returns a data.frame with the summary
#' @export
db_summary <- function(db) {
  if (!inherits(db, "RDatabase")) {
    stop("db must be an RDatabase object (created with vdjdb_open_file)")
  }

  cols <- db$to_columns()
  df <- as.data.frame(cols, stringsAsFactors = FALSE)

  cat("VDJdb Database Summary\n")
  cat("======================\n")
  cat("Total entries:", nrow(df), "\n\n")

  cat("By species:\n")
  species_counts <- sort(table(df$species), decreasing = TRUE)
  print(species_counts)

  cat("\nBy gene:\n")
  gene_counts <- sort(table(df$gene), decreasing = TRUE)
  print(gene_counts)

  cat("\nBy species and gene:\n")
  species_gene <- as.data.frame(table(df$species, df$gene))
  colnames(species_gene) <- c("species", "gene", "count")
  species_gene <- species_gene[species_gene$count > 0, ]
  species_gene <- species_gene[order(-species_gene$count), ]
  print(species_gene, row.names = FALSE)

  invisible(species_gene)
}
