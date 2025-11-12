#![allow(dead_code)]
use crate::error::{Result, VdjMatchError};
use crate::sequence::Clonotype;
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
        
        let mut entries = Vec::new();
        
        for result in reader.records() {
            let record = result?;
            
            // Parse record into DatabaseEntry
            // Note: This is simplified - real implementation should handle all VDJdb columns
        if record.len() >= 8 {
            let entry = DatabaseEntry {
                cdr3: record.get(1).unwrap_or("").to_string(),           // Column 1
                v_segment: record.get(7).unwrap_or("").to_string(),      // Column 7 (v.segm)
                j_segment: record.get(8).unwrap_or("").to_string(),      // Column 8 (j.segm)
                species: record.get(2).unwrap_or("").to_string(),        // Column 2
                gene: record.get(0).unwrap_or("").to_string(),           // Column 0
                mhc_class: Some(record.get(11).unwrap_or("").to_string()), // Column 11
                antigen_epitope: record.get(3).unwrap_or("").to_string(), // Column 3
                antigen_gene: Some(record.get(4).unwrap_or("").to_string()), // Column 4
                antigen_species: record.get(5).unwrap_or("").to_string(), // Column 5
                reference_id: record.get(12).map(|s| s.to_string()),     // Column 12
                method: None,  // Not in slim version
                meta: None,    // Not in slim version
                cdr3_fix: None, // Not in slim version
                vdjdb_score: record.get(13)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(0),
            };
                entries.push(entry);
            }
        }
        
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
