use crate::alignment::{align, matches_within_scope};
use crate::database::{Database, DatabaseEntry};
use crate::scoring::{compute_normalized_score, segment_match_score, simple_mismatch_score};
use crate::sequence::{Clonotype, SearchScope};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

/// A match between a query clonotype and a database entry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ClonotypeMatch {
    pub query_clonotype: Clonotype,
    pub db_entry: DatabaseEntry,
    pub score: f64,
    pub weight: f64,
    pub cdr3_alignment_score: f64,
    pub v_score: f64,
    pub j_score: f64,
    pub edit_distance: usize,
}

/// Configuration for matching
#[derive(Debug, Clone)]
pub struct MatchConfig {
    pub search_scope: SearchScope,
    pub match_v: bool,
    pub match_j: bool,
    pub use_vdjmatch_scoring: bool,
    pub scoring_mode: u8,
    pub exhaustive_search: u8,
    pub score_threshold: Option<f64>,
    pub max_hits_only: bool,
    pub top_n_hits: Option<usize>,
    pub weight_by_informativeness: bool,
}

impl Default for MatchConfig {
    fn default() -> Self {
        Self {
            search_scope: SearchScope::EXACT,
            match_v: false,
            match_j: false,
            use_vdjmatch_scoring: false,
            scoring_mode: 1,
            exhaustive_search: 1,
            score_threshold: None,
            max_hits_only: false,
            top_n_hits: None,
            weight_by_informativeness: false,
        }
    }
}

/// Match a clonotype against the database
pub fn match_clonotype(
    clonotype: &Clonotype,
    database: &Database,
    config: &MatchConfig,
) -> Vec<ClonotypeMatch> {
    let mut matches = Vec::new();
    
    for db_entry in &database.entries {
        // Check segment matches if required and if query has non-empty segments
        // Skip segment matching if query segment is empty (user wants CDR3-only matching)
        if config.match_v && !clonotype.v_segment.is_empty() {
            let v_match = Clonotype::normalize_segment(&clonotype.v_segment)
                == Clonotype::normalize_segment(&db_entry.v_segment);
            if !v_match {
                continue;
            }
        }

        if config.match_j && !clonotype.j_segment.is_empty() {
            let j_match = Clonotype::normalize_segment(&clonotype.j_segment)
                == Clonotype::normalize_segment(&db_entry.j_segment);
            if !j_match {
                continue;
            }
        }
        
        // Check CDR3 sequence match within scope
        let query_cdr3_str = &clonotype.cdr3_aa.sequence;
        let db_cdr3_str = &db_entry.cdr3;
        
        if !matches_within_scope(
            &clonotype.cdr3_aa,
            &crate::sequence::Cdr3Sequence::new(db_cdr3_str.clone()),
            &config.search_scope,
        ) {
            continue;
        }
        
        // Perform alignment
        let alignment = align(query_cdr3_str, db_cdr3_str);
        
        // Compute scores
        let cdr3_score = if config.use_vdjmatch_scoring {
            if config.scoring_mode == 1 {
                compute_normalized_score(&alignment)
            } else {
                simple_mismatch_score(&alignment)
            }
        } else {
            simple_mismatch_score(&alignment)
        };
        
        let v_score = segment_match_score(&clonotype.v_segment, &db_entry.v_segment, true);
        let j_score = segment_match_score(&clonotype.j_segment, &db_entry.j_segment, true);
        
        // Aggregate score
        let total_score = if config.use_vdjmatch_scoring {
            // VDJMATCH scoring: weighted combination
            0.5 * cdr3_score + 0.25 * v_score + 0.25 * j_score
        } else {
            cdr3_score
        };
        
        // Apply score threshold
        if let Some(threshold) = config.score_threshold {
            if total_score < threshold {
                continue;
            }
        }
        
        let matched = ClonotypeMatch {
            query_clonotype: clonotype.clone(),
            db_entry: db_entry.clone(),
            score: total_score,
            weight: 1.0, // Will be computed later if needed
            cdr3_alignment_score: cdr3_score,
            v_score,
            j_score,
            edit_distance: alignment.edit_distance,
        };
        
        matches.push(matched);
    }
    
    // Apply hit filtering
    if config.max_hits_only && !matches.is_empty() {
        let max_score = matches.iter().map(|m| m.score).fold(f64::NEG_INFINITY, f64::max);
        matches.retain(|m| (m.score - max_score).abs() < 1e-9);
    }
    
    if let Some(top_n) = config.top_n_hits {
        matches.sort_by(|a, b| b.score.partial_cmp(&a.score).unwrap());
        matches.truncate(top_n);
    }
    
    // Compute weights if requested
    if config.weight_by_informativeness {
        compute_informativeness_weights(&mut matches, database);
    }
    
    matches
}

/// Match multiple clonotypes in parallel
pub fn match_clonotypes_parallel(
    clonotypes: &[Clonotype],
    database: &Database,
    config: &MatchConfig,
) -> Vec<Vec<ClonotypeMatch>> {
    clonotypes
        .par_iter()
        .map(|clonotype| match_clonotype(clonotype, database, config))
        .collect()
}

/// Compute informativeness weights for matches
/// Weight = -log10(P(match by chance))
fn compute_informativeness_weights(matches: &mut [ClonotypeMatch], database: &Database) {
    // Group by epitope
    let mut epitope_cdr3_counts: std::collections::HashMap<String, usize> =
        std::collections::HashMap::new();
    
    for entry in &database.entries {
        *epitope_cdr3_counts
            .entry(entry.antigen_epitope.clone())
            .or_insert(0) += 1;
    }
    
    for m in matches.iter_mut() {
        let epitope = &m.db_entry.antigen_epitope;
        let count = epitope_cdr3_counts.get(epitope).copied().unwrap_or(1);
        
        // Simple informativeness: inverse of frequency
        // More sophisticated: based on probability of random match
        let total_entries = database.entries.len() as f64;
        let prob = (count as f64 + 1.0) / (total_entries + 1.0);
        m.weight = -prob.log10();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_match_clonotype() {
        let clonotype = Clonotype::new(
            "CASSLGQAYEQYF".to_string(),
            "TRBV12-3".to_string(),
            "TRBJ2-7".to_string(),
            100,
            0.01,
        );
        
        let db_entry = DatabaseEntry {
            cdr3: "CASSLGQAYEQYF".to_string(),
            v_segment: "TRBV12-3".to_string(),
            j_segment: "TRBJ2-7".to_string(),
            species: "HomoSapiens".to_string(),
            gene: "TRB".to_string(),
            mhc_class: Some("MHCI".to_string()),
            antigen_epitope: "GLCTLVAML".to_string(),
            antigen_gene: Some("BMLF1".to_string()),
            antigen_species: "EBV".to_string(),
            reference_id: Some("PMID:12345".to_string()),
            method: None,
            meta: None,
            cdr3_fix: None,
            vdjdb_score: 3,
        };
        
        let database = Database {
            entries: vec![db_entry],
            metadata: crate::database::DatabaseMetadata {
                columns: vec![],
                version: None,
            },
        };
        
        let config = MatchConfig::default();
        
        let matches = match_clonotype(&clonotype, &database, &config);
        
        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].score, 1.0);
    }
}
