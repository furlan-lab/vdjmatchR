#!/usr/bin/env Rscript
# Test script for db_to_table functionality

library(vdjmatchR)

cat("=== Testing db_to_table() Functionality ===\n\n")

# Test 1: Load database and get summary
cat("Test 1: Database Summary\n")
cat("------------------------\n")
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))
db_summary(db)
cat("\n")

# Test 2: Convert to data.table
cat("\nTest 2: Convert to data.table\n")
cat("------------------------------\n")

if (requireNamespace("data.table", quietly = TRUE)) {
  dt <- db_to_table(db)
  cat("Converted database to data.table\n")
  cat("Class:", class(dt), "\n")
  cat("Dimensions:", nrow(dt), "rows x", ncol(dt), "columns\n")
  cat("Columns:", paste(colnames(dt), collapse = ", "), "\n\n")

  # Show first few rows
  cat("First 5 rows:\n")
  print(head(dt, 5))
  cat("\n")

  # Test 3: Use data.table operations
  cat("\nTest 3: data.table Operations\n")
  cat("------------------------------\n")

  library(data.table)

  # Top 10 most common epitopes
  cat("Top 10 most common epitopes:\n")
  top_epitopes <- dt[, .N, by = antigen_epitope][order(-N)][1:10]
  print(top_epitopes)
  cat("\n")

  # TRB only
  cat("Human TRB entries:\n")
  dt_trb <- dt[species == "HomoSapiens" & gene == "TRB"]
  cat("  Total:", nrow(dt_trb), "entries\n")

  # Top TRB epitopes
  cat("\nTop 10 TRB epitopes:\n")
  top_trb <- dt_trb[, .N, by = antigen_epitope][order(-N)][1:10]
  print(top_trb)
  cat("\n")

  # V gene usage
  cat("TRB V gene usage (top 10):\n")
  v_usage <- dt_trb[, .N, by = v_segment][order(-N)][1:10]
  print(v_usage)
  cat("\n")

  # Test 4: Filter and convert
  cat("\nTest 4: Filter then Convert\n")
  cat("---------------------------\n")
  db_filtered <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)
  dt_filtered <- db_to_table(db_filtered)

  cat("Filtered to HomoSapiens TRB with vdjdb_score >= 2\n")
  cat("Entries:", nrow(dt_filtered), "\n")
  cat("Score distribution:\n")
  print(table(dt_filtered$vdjdb_score))
  cat("\n")

  # Test 5: Find specific epitopes
  cat("\nTest 5: Search for Specific Epitopes\n")
  cat("-------------------------------------\n")

  # GILGFVFTL is a common influenza epitope
  flu_hits <- dt_trb[antigen_epitope == "GILGFVFTL"]
  cat("Entries for GILGFVFTL (Influenza M1):", nrow(flu_hits), "\n")
  if (nrow(flu_hits) > 0) {
    cat("Unique CDR3 sequences:", uniqueN(flu_hits$cdr3), "\n")
    cat("Sample CDR3s:\n")
    print(head(unique(flu_hits$cdr3), 5))
  }
  cat("\n")

  # Test 6: Custom filtering
  cat("\nTest 6: Complex Filtering\n")
  cat("-------------------------\n")

  # Find epitopes with many different TCRs (public epitopes)
  epitope_diversity <- dt_trb[, .(
    n_entries = .N,
    n_unique_cdr3 = uniqueN(cdr3),
    n_unique_v = uniqueN(v_segment)
  ), by = antigen_epitope][order(-n_unique_cdr3)]

  cat("Top 5 epitopes by CDR3 diversity:\n")
  print(head(epitope_diversity, 5))
  cat("\n")

} else {
  cat("data.table not installed - testing with data.frame\n")
  df <- db_to_df(db)
  cat("Converted database to data.frame\n")
  cat("Class:", class(df), "\n")
  cat("Dimensions:", nrow(df), "rows x", ncol(df), "columns\n")
  cat("\nFirst 5 rows:\n")
  print(head(df, 5))
  cat("\nInstall data.table for advanced features:\n")
  cat("  install.packages('data.table')\n\n")
}

cat("=== All Tests Completed ===\n")
cat("\nKey Features:\n")
cat("✓ db_to_table() - Convert database to data.table/data.frame\n")
cat("✓ db_to_df() - Convert to data.frame explicitly\n")
cat("✓ db_summary() - Quick database overview\n")
cat("✓ Fast column-based extraction from Rust\n")
cat("✓ Full data.table compatibility for filtering/aggregation\n")
