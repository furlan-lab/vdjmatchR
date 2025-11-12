// Reuse core modules ported from vdjmatch-rs
pub mod alignment;
pub mod database;
pub mod error;
pub mod filtering;
pub mod matching;
pub mod scoring;
pub mod sequence;
pub mod utils;

use extendr_api::prelude::*;

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
        epitope.push(m.db_entry.antigen_epitope);
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
        epitope = epitope,
        score = score,
        cdr3_score = cdr3_score,
        v_score = v_score,
        j_score = j_score,
        edit_distance = edit_distance
    )
}

/// Batch match: vectors of cdr3/v/j; returns stacked results with query metadata.
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
    let mut score = Vec::new();
    let mut cdr3_score = Vec::new();
    let mut v_score = Vec::new();
    let mut j_score = Vec::new();
    let mut edit_distance = Vec::new();

    for (i, (cdr3i, (vi, ji))) in cdr3.into_iter().zip(v_segment.into_iter().zip(j_segment.into_iter())).enumerate() {
        let clonotype = sequence::Clonotype::new(cdr3i.clone(), vi.clone(), ji.clone(), 1, 0.0);
        let mut config = matching::MatchConfig::default();
        config.search_scope = search_scope;
        config.match_v = !vi.is_empty();
        config.match_j = !ji.is_empty();
        if top_n > 0 { config.top_n_hits = Some(top_n as usize); }

        let hits = matching::match_clonotype(&clonotype, &db.inner, &config);
        for m in hits.into_iter() {
            all_query_index.push((i as i32) + 1); // 1-based index for R
            all_query_cdr3.push(clonotype.cdr3_aa.sequence.clone());
            all_query_v.push(clonotype.v_segment.clone());
            all_query_j.push(clonotype.j_segment.clone());

            cdr3_db.push(m.db_entry.cdr3);
            v_db.push(m.db_entry.v_segment);
            j_db.push(m.db_entry.j_segment);
            species.push(m.db_entry.species);
            gene.push(m.db_entry.gene);
            epitope.push(m.db_entry.antigen_epitope);
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
        epitope = epitope,
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

// Register exported functions/types with R.
extendr_module! {
    mod vdjmatchR;
    impl RDatabase;
    fn match_tcr;
    fn match_tcr_many;
    fn vdjdb_ensure;
    fn vdjdb_update;
    fn vdjdb_ensure_into;
    fn vdjdb_update_into;
}
