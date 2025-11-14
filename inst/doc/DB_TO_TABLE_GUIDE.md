# Database to data.table/data.frame Guide

## Overview

New functions to convert the VDJdb database from Rust to R data structures for easy manipulation with standard R tools.

## Functions

### `db_to_table(db)`
Convert database to data.table (or data.frame if data.table not installed).

**Returns:** data.table with columns:
- `gene` - TRA or TRB
- `cdr3` - CDR3 amino acid sequence
- `v_segment` - V gene (e.g., "TRBV19*01")
- `j_segment` - J gene (e.g., "TRBJ2-1*01")
- `species` - Species (e.g., "HomoSapiens")
- `antigen_epitope` - Peptide sequence
- `antigen_gene` - Source gene of epitope
- `antigen_species` - Source species
- `mhc_class` - MHC restriction (MHCI/MHCII)
- `reference_id` - PubMed ID or reference
- `vdjdb_score` - Confidence score (0-3)

### `db_to_df(db)`
Convert database to data.frame explicitly.

### `db_summary(db)`
Print quick summary of database contents by species and gene.

## Basic Usage

```r
library(vdjmatchR)

# Load database
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))

# Quick summary
db_summary(db)

# Convert to data.table
dt <- db_to_table(db)

# Or data.frame
df <- db_to_df(db)
```

## Example Workflows

### 1. Explore Database Contents

```r
library(vdjmatchR)
library(data.table)

db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))
dt <- db_to_table(db)

# Most common epitopes
dt[, .N, by = antigen_epitope][order(-N)][1:10]

# Species distribution
dt[, .N, by = species]

# Gene distribution
dt[, .N, by = gene]

# Score quality
table(dt$vdjdb_score)
```

### 2. Filter Before Converting

```r
# Filter in Rust (faster for large filters)
db_trb <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)
dt_trb <- db_to_table(db_trb)

# Or filter in R
dt <- db_to_table(db)
dt_trb <- dt[species == "HomoSapiens" & gene == "TRB" & vdjdb_score >= 2]
```

### 3. Find TCRs for Specific Epitopes

```r
# Get all TCRs for a specific epitope
flu_tcrs <- dt[antigen_epitope == "GILGFVFTL"]  # Influenza M1

# Count unique CDR3s
flu_tcrs[, uniqueN(cdr3)]

# Get V gene usage for this epitope
flu_tcrs[, .N, by = v_segment][order(-N)]

# Find public epitopes (many different TCRs)
epitope_diversity <- dt[, .(
  n_tcrs = .N,
  n_unique_cdr3 = uniqueN(cdr3),
  avg_score = mean(vdjdb_score)
), by = antigen_epitope][order(-n_unique_cdr3)]
```

### 4. Create Custom Databases

```r
# Extract only high-confidence viral epitopes
viral_db <- dt[
  antigen_species %in% c("CMV", "EBV", "InfluenzaA", "SARS-CoV-2") &
  vdjdb_score >= 2
]

# Save as custom database
fwrite(viral_db, "viral_epitopes_high_conf.tsv", sep = "\t")

# Or filter by MHC restriction
mhci_db <- dt[mhc_class == "MHCI"]
```

### 5. TCR Repertoire Analysis

```r
# Find which epitopes your TCRs match
my_cdr3s <- c("CASSTHTNSYNEQFF", "CASSLGQAYEQYF", "CASSPGQGTTDTQYF")

# Exact matches
dt[cdr3 %in% my_cdr3s, .(cdr3, antigen_epitope, antigen_species)]

# Fuzzy matches (use vdjmatchR for this)
# But can explore similar sequences in data.table
dt[cdr3 %like% "CASS.*EQFF", .(cdr3, antigen_epitope)][1:10]
```

### 6. Statistical Analysis

```r
# CDR3 length distribution by gene
dt[, .(
  mean_len = mean(nchar(cdr3)),
  median_len = median(nchar(cdr3)),
  sd_len = sd(nchar(cdr3))
), by = gene]

# V gene usage frequency
v_usage <- dt[, .N, by = .(gene, v_segment)][order(-N)]

# Epitope size distribution
dt[, .(epitope_size = .N), by = antigen_epitope][, .(
  mean = mean(epitope_size),
  median = median(epitope_size),
  max = max(epitope_size)
)]

# Species with most epitopes
dt[, uniqueN(antigen_epitope), by = antigen_species][order(-V1)]
```

### 7. Join with Your Data

```r
library(data.table)

# Your experimental TCRs
my_tcrs <- data.table(
  sample_id = c("S1", "S1", "S2"),
  cell_id = c("cell1", "cell2", "cell3"),
  cdr3 = c("CASSTHTNSYNEQFF", "CASSLGQAYEQYF", "CASSPGQGTTDTQYF"),
  v_gene = c("TRBV19", "TRBV12-3", "TRBV5-4")
)

# Get database as data.table
dt <- db_to_table(db)

# Join to find exact matches
setkey(my_tcrs, cdr3)
setkey(dt, cdr3)
matched <- dt[my_tcrs, nomatch = 0]

# Or use vdjmatchR for fuzzy matching (recommended)
```

### 8. Create Filtered Databases for Matching

```r
# Create a tumor-specific database
tumor_epitopes <- c("MART1", "MAGE-A3", "NY-ESO-1", "gp100")

dt_tumor <- dt[antigen_gene %in% tumor_epitopes]

# Save filtered entries as new file
fwrite(dt_tumor, "tumor_epitopes.tsv", sep = "\t")

# Load back as database (would need proper format)
# Or use dt_tumor directly for lookups
```

## Performance Tips

1. **Filter in Rust first** (faster for large filters):
   ```r
   db_filtered <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)
   dt <- db_to_table(db_filtered)
   ```

2. **Use data.table** for better performance:
   ```r
   # Fast operations
   dt[, .N, by = gene]  # Much faster than aggregate()
   ```

3. **Convert once, reuse**:
   ```r
   dt <- db_to_table(db)  # Convert once
   # Multiple operations on dt...
   ```

4. **For very large databases**, filter before converting:
   ```r
   # Don't do this with full DB
   dt_full <- db_to_table(db)  # 214k rows
   dt_subset <- dt_full[species == "HomoSapiens"]  # Slow

   # Do this instead
   db_human <- filter_db(db, species = "HomoSapiens", gene = NULL, min_vdjdb_score = 0)
   dt_subset <- db_to_table(db_human)  # Fast
   ```

## Combining with TCR Matching

```r
library(vdjmatchR)
library(data.table)

# 1. Load database
db <- vdjdb_open_file(vdjdb_packaged_path(use_fat_db = TRUE))

# 2. Explore epitopes of interest
dt <- db_to_table(db)
viral_epitopes <- dt[
  antigen_species %in% c("CMV", "EBV") &
  vdjdb_score >= 2,
  .(antigen_epitope, antigen_species, .N),
  by = antigen_epitope
][order(-N)]

# 3. Filter database for matching
db_viral <- filter_db(db, species = "HomoSapiens", gene = "TRB", min_vdjdb_score = 2)

# 4. Match your TCRs
my_tcrs <- data.table(
  cdr3 = c("CASSTHTNSYNEQFF", "CASSLGQAYEQYF"),
  v = c("TRBV19", "TRBV12-3"),
  j = c("TRBJ2-1", "TRBJ2-7")
)

hits <- match_tcr_many_df(
  db_viral,
  my_tcrs$cdr3,
  my_tcrs$v,
  my_tcrs$j,
  scope = "3,1,2,3",  # Fuzzy matching
  top_n = 5
)

# 5. Annotate results
setDT(hits)
hits[, .(query_cdr3, epitope, epitope_species = species, score)][order(-score)]
```

## Summary

The `db_to_table()` family of functions provides:

✅ **Fast extraction** - Rust-based column extraction
✅ **Full R compatibility** - Works with data.table, dplyr, base R
✅ **Flexible filtering** - Filter before or after conversion
✅ **Easy exploration** - Standard data manipulation tools
✅ **Integration** - Combine with TCR matching workflows

Use these functions when you need to:
- Explore database contents
- Create custom filtered databases
- Join with your own data
- Perform statistical analyses
- Find specific epitopes or TCR patterns
