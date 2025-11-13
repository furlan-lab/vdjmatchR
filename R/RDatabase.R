# #' RDatabase: Handle to a loaded VDJdb
# #'
# #' @description
# #' `RDatabase` is an object provided by the Rust backend (via extendr) that
# #' holds an in-memory VDJdb table. Use `vdjdb_open_file()` (or helpers like
# #' `vdjdb_open()`) to construct a handle, `vdjdb_len()` to inspect how many
# #' rows were loaded, and `filter_db()` / `filter_db_by_epitope_size()` to derive
# #' filtered handles for downstream calls to \link[=match_tcr_df]{match_tcr_df}
# #' or \link[=match_tcr_many_df]{match_tcr_many_df}.
# #'
# #' See also: \link[=vdjdb_download]{vdjdb_download} for updating packaged
# #' databases.
# #'
# #' @name RDatabase
# #' @docType class
# #' @keywords internal
# NULL
