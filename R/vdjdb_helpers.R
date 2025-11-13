#' Open the packaged VDJdb via the Rust backend
#'
#' @inheritParams vdjdb_packaged_path
#' @inheritParams vdjdb_open
#' @inherit vdjdb_open_file return
#' @export
vdjdb_open_packaged <- function(use_fat_db = FALSE) {
  path <- vdjdb_packaged_path(use_fat_db)
  if (!nzchar(path) || !file.exists(path)) {
    stop("Packaged VDJdb file not found; set a user DB via vdjdb_set_user_db().")
  }
  vdjdb_open_file(path)
}

#' Open the preferred VDJdb (user path if set, else packaged)
#'
#' @inheritParams vdjdb_path
#' @return An in-memory VDJdb handle created by the Rust backend.
#' @export
vdjdb_open <- function(use_fat_db = FALSE) {
  path <- vdjdb_path(use_fat_db)
  vdjdb_open_file(path)
}
