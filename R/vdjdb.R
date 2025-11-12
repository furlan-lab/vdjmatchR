#' Locate packaged VDJdb file bundled with vdjmatchR
#'
#' @param use_fat_db logical; TRUE for full db if bundled, FALSE for slim (default)
#' @return string file path (empty string if not found)
#' @export
vdjdb_packaged_path <- function(use_fat_db = FALSE) {
  # Prefer an updated uncompressed file if present (downloaded into extdata)
  fname <- if (isTRUE(use_fat_db)) "vdjdb.txt" else "vdjdb.slim.txt"
  path_txt <- system.file("extdata", fname, package = "vdjmatchR", mustWork = FALSE)
  if (!is.null(path_txt) && nzchar(path_txt) && file.exists(path_txt)) return(path_txt)

  # Fall back to gzipped asset bundled with the package
  gz <- if (isTRUE(use_fat_db)) "vdjdb.txt.gz" else "vdjdb.slim.txt.gz"
  path_gz <- system.file("extdata", gz, package = "vdjmatchR", mustWork = FALSE)
  if (!is.null(path_gz) && nzchar(path_gz)) return(path_gz)
  ""
}

#' Preferred VDJdb path (user-specified if set, else packaged)
#'
#' Returns a user-specified path if configured via \link[=vdjdb_set_user_db]{vdjdb_set_user_db},
#' otherwise returns the packaged file path.
#'
#' @param use_fat_db logical; TRUE for full db, FALSE for slim (default)
#' @return string file path
#' @export
vdjdb_path <- function(use_fat_db = FALSE) {
  # Prefer user-specified DB path set via options
  opt_name <- if (isTRUE(use_fat_db)) "vdjmatchR.vdjdb_path_fat" else "vdjmatchR.vdjdb_path_slim"
  user_path <- getOption(opt_name, default = "")
  if (nzchar(user_path)) {
    if (!file.exists(user_path)) stop(sprintf("User-specified VDJdb file not found: %s", user_path))
    return(user_path)
  }

  # Fallback to packaged file
  pkg <- vdjdb_packaged_path(use_fat_db)
  if (nzchar(pkg)) return(pkg)

  stop("No VDJdb file found: provide a file via vdjdb_set_user_db() or use the packaged one.")
}

#' Update (download) the VDJdb file and return its path into extdata
#'
#' @param use_fat_db logical; TRUE for full db, FALSE for slim (default)
#' @return string file path to the updated file
#' @export
vdjdb_update_latest <- function(use_fat_db = FALSE) {
  # Save the updated file into the installed package's extdata directory
  data_dir <- system.file("extdata", package = "vdjmatchR")
  if (!nzchar(data_dir)) stop("Could not resolve package extdata directory.")
  # Check writability of extdata directory
  probe <- file.path(data_dir, ".vdjmatchR_write_test")
  can_write <- isTRUE(tryCatch({ file.create(probe); file.remove(probe); TRUE }, error = function(e) FALSE))
  if (!can_write) {
    stop(sprintf("Package extdata directory is not writable: %s. Set a user-supplied database with vdjdb_set_user_db().", data_dir))
  }
  # Trigger ensure/download of the requested variant into extdata
  p <- vdjdb_ensure_into(data_dir, isTRUE(use_fat_db))
  return(p)
}

#' Back-compat alias for updating/downloading VDJdb (one variant)
#' @inheritParams vdjdb_update_latest
#' @inherit vdjdb_update_latest return
#' @export
vdjdb_download <- function(use_fat_db = FALSE) {
  vdjdb_update_latest(use_fat_db)
}

#' Set a user-specified VDJdb path to use by default
#'
#' Stores the path in R options for the session. Call with `use_fat_db = TRUE`
#' to set the full database, or FALSE for the slim one.
#'
#' @param path Path to a VDJdb TSV/TSV.GZ file
#' @param use_fat_db logical; TRUE for full db, FALSE for slim (default)
#' @export
vdjdb_set_user_db <- function(path, use_fat_db = FALSE) {
  if (!file.exists(path)) stop(sprintf("File does not exist: %s", path))
  opt_name <- if (isTRUE(use_fat_db)) "vdjmatchR.vdjdb_path_fat" else "vdjmatchR.vdjdb_path_slim"
  options(structure(list(path), .Names = opt_name))
  invisible(normalizePath(path))
}

#' Update/download both slim and fat VDJdb files (disabled)
#' @export
vdjdb_update_all <- function() {
  data_dir <- system.file("extdata", package = "vdjmatchR")
  if (!nzchar(data_dir)) stop("Could not resolve package extdata directory.")
  # Check writability of extdata directory
  probe <- file.path(data_dir, ".vdjmatchR_write_test")
  can_write <- isTRUE(tryCatch({ file.create(probe); file.remove(probe); TRUE }, error = function(e) FALSE))
  if (!can_write) {
    stop(sprintf("Package extdata directory is not writable: %s. Set a user-supplied database with vdjdb_set_user_db().", data_dir))
  }
  vdjdb_update_into(data_dir)
  invisible(NULL)
}
