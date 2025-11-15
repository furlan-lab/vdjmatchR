#![allow(non_snake_case)]

// Reuse core modules ported from vdjmatch-rs
pub mod alignment;
pub mod database;
pub mod error;
pub mod filtering;
pub mod matching;
pub mod scoring;
pub mod sequence;
pub mod tcrdist;
pub mod utils;

use extendr_api::prelude::*;
use std::path::Path;

#[extendr]
pub struct RDatabase {
    inner: database::Database,
}

#[extendr]
impl RDatabase {
    pub fn new_from_file(path: &str) -> Result<Self> {
        match database::Database::load_from_file(path) {
            Ok(db) => Ok(Self { inner: db }),
            Err(e) => Err(extendr_api::error::Error::Other(e.to_string())),
        }
    }

    pub fn new_from_vdjdb(_use_fat_db: bool) -> Result<Self> {
        Err(extendr_api::error::Error::Other(
            "new_from_vdjdb is disabled. Provide an explicit file path via new_from_file().".into(),
        ))
    }

    pub fn len(&self) -> i32 {
        self.inner.len() as i32
    }

    /// Return a filtered copy of the database. Use NULL for no filter.
    pub fn filter(&self, species: Option<String>, gene: Option<String>, min_vdjdb_score: i32) -> Self {
        let filtered = self.inner.filter(
            species.as_deref(),
            gene.as_deref(),
            min_vdjdb_score as u8,
        );
        Self { inner: filtered }
    }

    /// Filter by minimum epitope size (unique CDR3s per epitope)
    pub fn filter_by_epitope_size(&self, min_size: i32) -> Self {
        let filtered = self.inner.filter_by_epitope_size(min_size as usize);
        Self { inner: filtered }
    }

    /// Convert database to column vectors for R data.frame/data.table
    pub fn to_columns(&self) -> List {
        let n = self.inner.entries.len();

        let mut gene = Vec::with_capacity(n);
        let mut cdr3 = Vec::with_capacity(n);
        let mut v_segment = Vec::with_capacity(n);
        let mut j_segment = Vec::with_capacity(n);
        let mut species = Vec::with_capacity(n);
        let mut antigen_epitope = Vec::with_capacity(n);
        let mut antigen_gene = Vec::with_capacity(n);
        let mut antigen_species = Vec::with_capacity(n);
        let mut mhc_class = Vec::with_capacity(n);
        let mut reference_id = Vec::with_capacity(n);
        let mut vdjdb_score = Vec::with_capacity(n);

        for entry in &self.inner.entries {
            gene.push(entry.gene.clone());
            cdr3.push(entry.cdr3.clone());
            v_segment.push(entry.v_segment.clone());
            j_segment.push(entry.j_segment.clone());
            species.push(entry.species.clone());
            antigen_epitope.push(entry.antigen_epitope.clone());
            antigen_gene.push(entry.antigen_gene.clone().unwrap_or_default());
            antigen_species.push(entry.antigen_species.clone());
            mhc_class.push(entry.mhc_class.clone().unwrap_or_default());
            reference_id.push(entry.reference_id.clone().unwrap_or_default());
            vdjdb_score.push(entry.vdjdb_score as i32);
        }

        list!(
            gene = gene,
            cdr3 = cdr3,
            v_segment = v_segment,
            j_segment = j_segment,
            species = species,
            antigen_epitope = antigen_epitope,
            antigen_gene = antigen_gene,
            antigen_species = antigen_species,
            mhc_class = mhc_class,
            reference_id = reference_id,
            vdjdb_score = vdjdb_score
        )
    }
}

/// Open a VDJdb TSV/TSV.GZ via the Rust backend.
/// @export
#[extendr]
pub fn vdjdb_open_file(path: &str) -> Result<RDatabase> {
    if path.trim().is_empty() {
        return Err(extendr_api::error::Error::Other("path must be a non-empty string".into()));
    }
    if !Path::new(path).exists() {
        return Err(extendr_api::error::Error::Other(format!("VDJdb file not found: {path}")));
    }
    RDatabase::new_from_file(path)
}

/// Number of rows stored in the in-memory VDJdb handle.
/// @export
#[extendr]
pub fn vdjdb_len(db: &RDatabase) -> i32 {
    db.len()
}

/// Filter database entries by species, gene, and minimum VDJdb score.
/// @export
#[extendr]
pub fn filter_db(
    db: &RDatabase,
    species: Nullable<String>,
    gene: Nullable<String>,
    min_vdjdb_score: i32,
) -> RDatabase {
    let species_string = species.into_option().filter(|s| !s.trim().is_empty());
    let gene_string = gene.into_option().filter(|s| !s.trim().is_empty());
    db.filter(species_string, gene_string, min_vdjdb_score)
}

/// Filter by minimum epitope size (unique CDR3 per epitope).
/// @export
#[extendr]
pub fn filter_db_by_epitope_size(db: &RDatabase, min_size: i32) -> RDatabase {
    db.filter_by_epitope_size(min_size)
}

/// Match a single clonotype against the database.
/// Returns a list of columns (vector-of-equal-length) suitable for as.data.frame in R.
#[extendr]
pub fn match_tcr(
    db: &RDatabase,
    cdr3: &str,
    v_segment: &str,
    j_segment: &str,
    scope: &str,
    top_n: i32,
) -> List {
    let clonotype = sequence::Clonotype::new(
        cdr3.to_string(),
        v_segment.to_string(),
        j_segment.to_string(),
        1,
        0.0,
    );

    // Parse scope, default to exact on failure.
    let search_scope = sequence::SearchScope::parse(scope).unwrap_or(sequence::SearchScope::EXACT);

    let mut config = matching::MatchConfig::default();
    config.search_scope = search_scope;
    config.match_v = !v_segment.is_empty();
    config.match_j = !j_segment.is_empty();
    if top_n > 0 { config.top_n_hits = Some(top_n as usize); }

    let matches = matching::match_clonotype(&clonotype, &db.inner, &config);

    let n = matches.len();
    let mut cdr3_db = Vec::with_capacity(n);
    let mut v_db = Vec::with_capacity(n);
    let mut j_db = Vec::with_capacity(n);
    let mut species = Vec::with_capacity(n);
    let mut gene = Vec::with_capacity(n);
    let mut epitope = Vec::with_capacity(n);
    let mut antigen_gene = Vec::with_capacity(n);
    let mut antigen_species = Vec::with_capacity(n);
    let mut mhc_class = Vec::with_capacity(n);
    let mut reference_id = Vec::with_capacity(n);
    let mut vdjdb_score = Vec::with_capacity(n);
    let mut score = Vec::with_capacity(n);
    let mut cdr3_score = Vec::with_capacity(n);
    let mut v_score = Vec::with_capacity(n);
    let mut j_score = Vec::with_capacity(n);
    let mut edit_distance = Vec::with_capacity(n);

    for m in matches.into_iter() {
        cdr3_db.push(m.db_entry.cdr3);
        v_db.push(m.db_entry.v_segment);
        j_db.push(m.db_entry.j_segment);
        species.push(m.db_entry.species);
        gene.push(m.db_entry.gene);
        epitope.push(m.db_entry.antigen_epitope.clone());
        antigen_gene.push(m.db_entry.antigen_gene.unwrap_or_default());
        antigen_species.push(m.db_entry.antigen_species);
        mhc_class.push(m.db_entry.mhc_class.unwrap_or_default());
        reference_id.push(m.db_entry.reference_id.unwrap_or_default());
        vdjdb_score.push(m.db_entry.vdjdb_score as i32);
        score.push(m.score);
        cdr3_score.push(m.cdr3_alignment_score);
        v_score.push(m.v_score);
        j_score.push(m.j_score);
        edit_distance.push(m.edit_distance as i32);
    }

    list!(
        cdr3_db = cdr3_db,
        v_db = v_db,
        j_db = j_db,
        species = species,
        gene = gene,
        antigen_epitope = epitope,
        antigen_gene = antigen_gene,
        antigen_species = antigen_species,
        mhc_class = mhc_class,
        reference_id = reference_id,
        vdjdb_score = vdjdb_score,
        score = score,
        cdr3_score = cdr3_score,
        v_score = v_score,
        j_score = j_score,
        edit_distance = edit_distance
    )
}

/// Batch match: vectors of cdr3/v/j; returns stacked results with query metadata.
/// Uses parallel processing via Rayon for improved performance.
#[extendr]
pub fn match_tcr_many(
    db: &RDatabase,
    cdr3: Vec<String>,
    v_segment: Vec<String>,
    j_segment: Vec<String>,
    scope: &str,
    top_n: i32,
) -> Result<List> {
    if !(cdr3.len() == v_segment.len() && v_segment.len() == j_segment.len()) {
        return Err(extendr_api::error::Error::Other("cdr3, v_segment, j_segment must have equal length".into()));
    }

    let search_scope = sequence::SearchScope::parse(scope).unwrap_or(sequence::SearchScope::EXACT);

    // Build clonotypes for parallel matching
    let clonotypes: Vec<sequence::Clonotype> = cdr3
        .iter()
        .zip(v_segment.iter().zip(j_segment.iter()))
        .map(|(cdr3i, (vi, ji))| {
            sequence::Clonotype::new(cdr3i.clone(), vi.clone(), ji.clone(), 1, 0.0)
        })
        .collect();

    // Configure matching
    let mut config = matching::MatchConfig::default();
    config.search_scope = search_scope;
    config.match_v = true;  // Matching logic handles empty segments
    config.match_j = true;  // Matching logic handles empty segments
    if top_n > 0 { config.top_n_hits = Some(top_n as usize); }

    // Use parallel matching
    let all_matches = matching::match_clonotypes_parallel(&clonotypes, &db.inner, &config);

    // Flatten results
    let mut all_query_index: Vec<i32> = Vec::new();
    let mut all_query_cdr3: Vec<String> = Vec::new();
    let mut all_query_v: Vec<String> = Vec::new();
    let mut all_query_j: Vec<String> = Vec::new();

    let mut cdr3_db = Vec::new();
    let mut v_db = Vec::new();
    let mut j_db = Vec::new();
    let mut species = Vec::new();
    let mut gene = Vec::new();
    let mut epitope = Vec::new();
    let mut antigen_gene = Vec::new();
    let mut antigen_species = Vec::new();
    let mut mhc_class = Vec::new();
    let mut reference_id = Vec::new();
    let mut vdjdb_score = Vec::new();
    let mut score = Vec::new();
    let mut cdr3_score = Vec::new();
    let mut v_score = Vec::new();
    let mut j_score = Vec::new();
    let mut edit_distance = Vec::new();

    for (i, matches) in all_matches.into_iter().enumerate() {
        let clonotype = &clonotypes[i];
        for m in matches.into_iter() {
            all_query_index.push((i as i32) + 1); // 1-based index for R
            all_query_cdr3.push(clonotype.cdr3_aa.sequence.clone());
            all_query_v.push(clonotype.v_segment.clone());
            all_query_j.push(clonotype.j_segment.clone());

            cdr3_db.push(m.db_entry.cdr3);
            v_db.push(m.db_entry.v_segment);
            j_db.push(m.db_entry.j_segment);
            species.push(m.db_entry.species);
            gene.push(m.db_entry.gene);
            epitope.push(m.db_entry.antigen_epitope.clone());
            antigen_gene.push(m.db_entry.antigen_gene.unwrap_or_default());
            antigen_species.push(m.db_entry.antigen_species);
            mhc_class.push(m.db_entry.mhc_class.unwrap_or_default());
            reference_id.push(m.db_entry.reference_id.unwrap_or_default());
            vdjdb_score.push(m.db_entry.vdjdb_score as i32);
            score.push(m.score);
            cdr3_score.push(m.cdr3_alignment_score);
            v_score.push(m.v_score);
            j_score.push(m.j_score);
            edit_distance.push(m.edit_distance as i32);
        }
    }

    Ok(list!(
        query_index = all_query_index,
        query_cdr3 = all_query_cdr3,
        query_v = all_query_v,
        query_j = all_query_j,
        cdr3_db = cdr3_db,
        v_db = v_db,
        j_db = j_db,
        species = species,
        gene = gene,
        antigen_epitope = epitope,
        antigen_gene = antigen_gene,
        antigen_species = antigen_species,
        mhc_class = mhc_class,
        reference_id = reference_id,
        vdjdb_score = vdjdb_score,
        score = score,
        cdr3_score = cdr3_score,
        v_score = v_score,
        j_score = j_score,
        edit_distance = edit_distance
    ))
}

/// Ensure VDJdb exists locally and return the path.
#[extendr]
pub fn vdjdb_ensure(_use_fat_db: bool) -> Result<String> {
    Err(extendr_api::error::Error::Other(
        "Automatic download to home directory is disabled. Use vdjdb_ensure_into(dir, use_fat_db).".into(),
    ))
}

/// Download/update the VDJdb files (slim and fat).
#[extendr]
pub fn vdjdb_update() -> Result<()> {
    Err(extendr_api::error::Error::Other(
        "Automatic update to home directory is disabled. Use vdjdb_update_into(dir).".into(),
    ))
}

/// Ensure VDJdb exists in the specified directory and return the path.
#[extendr]
pub fn vdjdb_ensure_into(dir: &str, use_fat_db: bool) -> Result<String> {
    let mgr = database::DatabaseManager::new_with_dir(dir);
    match mgr.ensure_database_exists(use_fat_db) {
        Ok(path) => Ok(path.to_string_lossy().to_string()),
        Err(e) => Err(extendr_api::error::Error::Other(e.to_string())),
    }
}

/// Download/update the VDJdb files (slim and fat) into the specified directory.
#[extendr]
pub fn vdjdb_update_into(dir: &str) -> Result<()> {
    let mgr = database::DatabaseManager::new_with_dir(dir);
    match mgr.update_database() {
        Ok(()) => Ok(()),
        Err(e) => Err(extendr_api::error::Error::Other(e.to_string())),
    }
}

/// Calculate pairwise tcrdist distances between TCRs
/// Returns a distance matrix (as a vector in column-major order for R)
/// Pass empty strings for missing CDR sequences
/// @export
#[extendr]
pub fn calculate_tcrdist(
    cdr1_a: Vec<String>,
    cdr2_a: Vec<String>,
    cdr3_a: Vec<String>,
    cdr1_b: Vec<String>,
    cdr2_b: Vec<String>,
    cdr3_b: Vec<String>,
) -> Result<List> {
    let n = cdr3_a.len();

    // Validate input lengths
    if !(cdr1_a.len() == n && cdr2_a.len() == n &&
         cdr1_b.len() == n && cdr2_b.len() == n && cdr3_b.len() == n) {
        return Err(extendr_api::error::Error::Other(
            "All CDR vectors must have equal length".into()
        ));
    }

    // Helper to convert empty string to None
    let to_opt = |s: &str| if s.is_empty() { None } else { Some(s.to_string()) };

    // Build TCR objects
    let tcrs: Vec<tcrdist::TCR> = (0..n).map(|i| {
        tcrdist::TCR::new(
            to_opt(&cdr1_a[i]),
            to_opt(&cdr2_a[i]),
            to_opt(&cdr3_a[i]),
            to_opt(&cdr1_b[i]),
            to_opt(&cdr2_b[i]),
            to_opt(&cdr3_b[i]),
        )
    }).collect();

    // Calculate pairwise distances
    let mut distances = Vec::with_capacity(n * n);
    let mut i_indices = Vec::with_capacity(n * n);
    let mut j_indices = Vec::with_capacity(n * n);

    for i in 0..n {
        for j in 0..n {
            let dist = tcrdist::tcrdist(&tcrs[i], &tcrs[j]);
            i_indices.push((i + 1) as i32); // 1-based for R
            j_indices.push((j + 1) as i32);
            distances.push(dist);
        }
    }

    Ok(list!(
        i = i_indices,
        j = j_indices,
        distance = distances,
        n = n as i32
    ))
}

/// Calculate tcrdist between two single TCRs
/// Pass empty strings for missing CDR sequences
#[extendr]
pub fn tcrdist_single(
    cdr1_a_1: &str,
    cdr2_a_1: &str,
    cdr3_a_1: &str,
    cdr1_b_1: &str,
    cdr2_b_1: &str,
    cdr3_b_1: &str,
    cdr1_a_2: &str,
    cdr2_a_2: &str,
    cdr3_a_2: &str,
    cdr1_b_2: &str,
    cdr2_b_2: &str,
    cdr3_b_2: &str,
) -> f64 {
    let to_opt = |s: &str| if s.is_empty() { None } else { Some(s.to_string()) };

    let tcr1 = tcrdist::TCR::new(
        to_opt(cdr1_a_1),
        to_opt(cdr2_a_1),
        to_opt(cdr3_a_1),
        to_opt(cdr1_b_1),
        to_opt(cdr2_b_1),
        to_opt(cdr3_b_1),
    );

    let tcr2 = tcrdist::TCR::new(
        to_opt(cdr1_a_2),
        to_opt(cdr2_a_2),
        to_opt(cdr3_a_2),
        to_opt(cdr1_b_2),
        to_opt(cdr2_b_2),
        to_opt(cdr3_b_2),
    );

    tcrdist::tcrdist(&tcr1, &tcr2)
}

// Register exported functions/types with R.
extendr_module! {
    mod vdjmatchR;
    impl RDatabase;
    fn match_tcr;
    fn match_tcr_many;
    fn vdjdb_open_file;
    fn vdjdb_len;
    fn filter_db;
    fn filter_db_by_epitope_size;
    fn vdjdb_ensure;
    fn vdjdb_update;
    fn vdjdb_ensure_into;
    fn vdjdb_update_into;
    fn calculate_tcrdist;
    fn tcrdist_single;
}
