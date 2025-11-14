#![allow(dead_code)]
use crate::error::{Result, VdjMatchError};
// use crate::sequence::Clonotype;
use csv::ReaderBuilder;
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};

/// VDJdb database entry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DatabaseEntry {
    pub cdr3: String,
    pub v_segment: String,
    pub j_segment: String,
    pub species: String,
    pub gene: String,
    pub mhc_class: Option<String>,
    pub antigen_epitope: String,
    pub antigen_gene: Option<String>,
    pub antigen_species: String,
    pub reference_id: Option<String>,
    pub method: Option<String>,
    pub meta: Option<String>,
    pub cdr3_fix: Option<String>,
    pub vdjdb_score: u8,
}

impl DatabaseEntry {
    /// Check if this entry matches the given filters
    pub fn matches_species(&self, species: &str) -> bool {
        self.species.eq_ignore_ascii_case(species)
    }
    
    pub fn matches_gene(&self, gene: &str) -> bool {
        self.gene.eq_ignore_ascii_case(gene)
    }
    
    pub fn matches_vdjdb_score(&self, min_score: u8) -> bool {
        self.vdjdb_score >= min_score
    }
}

/// VDJdb database manager
pub struct Database {
    pub entries: Vec<DatabaseEntry>,
    pub metadata: DatabaseMetadata,
}

#[derive(Debug, Clone)]
pub struct DatabaseMetadata {
    pub columns: Vec<String>,
    pub version: Option<String>,
}

impl Database {
    /// Load database from file
    pub fn load_from_file<P: AsRef<Path>>(path: P) -> Result<Self> {
        let p = path.as_ref();
        let file = File::open(p)
            .map_err(|e| VdjMatchError::DatabaseNotFound(e.to_string()))?;

        // Support gzip-compressed TSVs by checking the file extension.
        let is_gz = p
            .extension()
            .and_then(|s| s.to_str())
            .map(|s| s.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);

        let reader: Box<dyn Read> = if is_gz {
            Box::new(GzDecoder::new(file))
        } else {
            Box::new(file)
        };

        let mut reader = ReaderBuilder::new()
            .delimiter(b'\t')
            .from_reader(BufReader::new(reader));
        
        let headers = reader.headers()?;
        let columns: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

        // Build column name -> index map for flexible column ordering
        let mut col_map = HashMap::new();
        for (i, col_name) in columns.iter().enumerate() {
            col_map.insert(col_name.as_str(), i);
        }

        // Required column names
        let gene_idx = col_map.get("gene").copied();
        let cdr3_idx = col_map.get("cdr3").copied();
        let species_idx = col_map.get("species").copied();
        let v_segm_idx = col_map.get("v.segm").copied();
        let j_segm_idx = col_map.get("j.segm").copied();
        let antigen_epitope_idx = col_map.get("antigen.epitope").copied();
        let antigen_gene_idx = col_map.get("antigen.gene").copied();
        let antigen_species_idx = col_map.get("antigen.species").copied();
        let mhc_class_idx = col_map.get("mhc.class").copied();
        let reference_id_idx = col_map.get("reference.id").copied();
        let vdjdb_score_idx = col_map.get("vdjdb.score").copied();
        let method_idx = col_map.get("method").copied();
        let meta_idx = col_map.get("meta").copied();
        let cdr3fix_idx = col_map.get("cdr3fix").copied();

        let mut entries = Vec::new();

        for result in reader.records() {
            let record = result?;

            // Parse record into DatabaseEntry using column names
            let entry = DatabaseEntry {
                gene: gene_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                cdr3: cdr3_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                v_segment: v_segm_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                j_segment: j_segm_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                species: species_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                antigen_epitope: antigen_epitope_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                antigen_gene: antigen_gene_idx.and_then(|i| record.get(i).map(|s| s.to_string())),
                antigen_species: antigen_species_idx.and_then(|i| record.get(i)).unwrap_or("").to_string(),
                mhc_class: mhc_class_idx.and_then(|i| record.get(i).map(|s| s.to_string())),
                reference_id: reference_id_idx.and_then(|i| record.get(i).map(|s| s.to_string())),
                method: method_idx.and_then(|i| record.get(i).map(|s| s.to_string())),
                meta: meta_idx.and_then(|i| record.get(i).map(|s| s.to_string())),
                cdr3_fix: cdr3fix_idx.and_then(|i| record.get(i).map(|s| s.to_string())),
                vdjdb_score: vdjdb_score_idx
                    .and_then(|i| record.get(i))
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0),
            };
            entries.push(entry);
        }

        // eprintln!("DEBUG: Loaded {} entries from database", entries.len());
        // if !entries.is_empty() {
        //     eprintln!("DEBUG: First 3 loaded entries:");
        //     for (i, entry) in entries.iter().take(3).enumerate() {
        //         eprintln!("  Entry {}: gene='{}' species='{}' cdr3='{}' v='{}' j='{}'",
        //                   i+1, entry.gene, entry.species, entry.cdr3, entry.v_segment, entry.j_segment);
        //     }
        // }

        Ok(Self {
            entries,
            metadata: DatabaseMetadata {
                columns,
                version: None,
            },
        })
    }
    
    /// Filter database entries by criteria
    pub fn filter(
        &self,
        species: Option<&str>,
        gene: Option<&str>,
        min_vdjdb_score: u8,
    ) -> Self {
        // eprintln!("DEBUG: Filtering {} entries", self.entries.len());
        // eprintln!("DEBUG: species filter={:?}, gene filter={:?}", species, gene);
        // if let Some(first) = self.entries.first() {
        //     eprintln!("DEBUG: First entry: gene='{}' species='{}'", first.gene, first.species);
        // }
        let filtered_entries: Vec<DatabaseEntry> = self
            .entries
            .iter()
            .filter(|entry| {
                if let Some(s) = species {
                    if !entry.matches_species(s) {
                        return false;
                    }
                }
                if let Some(g) = gene {
                    if !entry.matches_gene(g) {
                        return false;
                    }
                }
                if !entry.matches_vdjdb_score(min_vdjdb_score) {
                    return false;
                }
                true
            })
            .cloned()
            .collect();

        // eprintln!("DEBUG: After filtering: {} entries", filtered_entries.len());
        // eprintln!("DEBUG: First 3 filtered entries:");
        // for (i, entry) in filtered_entries.iter().take(3).enumerate() {
        //     eprintln!("  Entry {}: gene='{}' species='{}' cdr3='{}'",
        //               i+1, entry.gene, entry.species, entry.cdr3);
        // }

        Self {
            entries: filtered_entries,
            metadata: self.metadata.clone(),
        }

    }
    
    /// Filter by epitope size (minimum number of unique CDR3 per epitope)
    pub fn filter_by_epitope_size(&self, min_size: usize) -> Self {
        let mut epitope_counts: HashMap<String, usize> = HashMap::new();
        
        // Count unique CDR3s per epitope
        for entry in &self.entries {
            *epitope_counts.entry(entry.antigen_epitope.clone()).or_insert(0) += 1;
        }
        
        let filtered_entries: Vec<DatabaseEntry> = self
            .entries
            .iter()
            .filter(|entry| {
                epitope_counts
                    .get(&entry.antigen_epitope)
                    .map(|&count| count >= min_size)
                    .unwrap_or(false)
            })
            .cloned()
            .collect();

        // eprintln!("DEBUG: First 3 entries:");
        // for (i, entry) in filtered_entries.iter().take(3).enumerate() {
        //     eprintln!("  Entry {}: gene='{}' species='{}' cdr3='{}'",
        //               i+1, entry.gene, entry.species, entry.cdr3);
        // }

        Self {
            entries: filtered_entries,
            metadata: self.metadata.clone(),
        }
    }
    
    pub fn len(&self) -> usize {
        self.entries.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Database downloader and manager
pub struct DatabaseManager {
    home_dir: PathBuf,
}

impl DatabaseManager {
    pub fn new() -> Self {
        let home = dirs::home_dir().unwrap_or_else(|| PathBuf::from("."));
        let vdjdb_home = home.join(".vdjmatch");
        
        Self { home_dir: vdjdb_home }
    }
    
    /// Create a manager that stores files in the specified directory
    pub fn new_with_dir<P: AsRef<Path>>(dir: P) -> Self {
        Self { home_dir: dir.as_ref().to_path_buf() }
    }
    
    pub fn ensure_database_exists(&self, use_fat_db: bool) -> Result<PathBuf> {
        std::fs::create_dir_all(&self.home_dir)?;
        
        let db_file = if use_fat_db {
            self.home_dir.join("vdjdb.txt")
        } else {
            self.home_dir.join("vdjdb.slim.txt")
        };
        
        if !db_file.exists() {
            eprintln!("Database not found. Downloading...");
            self.download_database(use_fat_db)?;
        }
        
        Ok(db_file)
    }
    
    fn download_database(&self, use_fat_db: bool) -> Result<()> {
        let url = if use_fat_db {
            "https://github.com/antigenomics/vdjdb-db/releases/latest/download/vdjdb.txt"
        } else {
            "https://github.com/antigenomics/vdjdb-db/releases/latest/download/vdjdb.slim.txt"
        };
        
        eprintln!("Downloading from: {}", url);
        
        let response = reqwest::blocking::get(url)?;
        let content = response.bytes()?;
        
        let db_file = if use_fat_db {
            self.home_dir.join("vdjdb.txt")
        } else {
            self.home_dir.join("vdjdb.slim.txt")
        };
        
        let mut file = File::create(&db_file)?;
        file.write_all(&content)?;
        
        eprintln!("Database downloaded successfully");
        
        Ok(())
    }
    
    pub fn update_database(&self) -> Result<()> {
        eprintln!("Updating VDJdb database...");
        
        // Ensure directory exists
        std::fs::create_dir_all(&self.home_dir)?;
        
        // Download both versions
        self.download_database(true)?;
        self.download_database(false)?;
        
        eprintln!("Database updated successfully");
        
        Ok(())
    }
}

impl Default for DatabaseManager {
    fn default() -> Self {
        Self::new()
    }
}
