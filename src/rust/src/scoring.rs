use crate::alignment::{Alignment, EditOp};
use std::collections::HashMap;

lazy_static::lazy_static! {
    /// BLOSUM62-based substitution matrix for amino acids
    static ref BLOSUM62: HashMap<(u8, u8), i32> = {
        let mut m = HashMap::new();
        
        // Simplified BLOSUM62 matrix for key amino acids in TCR CDR3
        let acids = b"ARNDCQEGHILKMFPSTWYV";
        
        // Diagonal (matches) - positive scores
        for &aa in acids {
            m.insert((aa, aa), 4);
        }
        
        // Similar amino acids - small positive/negative scores
        let similar = vec![
            (b'A', b'S', 1), (b'A', b'T', 0),
            (b'R', b'K', 2), (b'R', b'Q', 1),
            (b'N', b'D', 1), (b'N', b'S', 1),
            (b'D', b'E', 2),
            (b'C', b'S', -1),
            (b'Q', b'E', 0), (b'Q', b'K', 1),
            (b'E', b'K', 1),
            (b'I', b'L', 2), (b'I', b'M', 1), (b'I', b'V', 3),
            (b'L', b'M', 2), (b'L', b'V', 1),
            (b'M', b'V', 1),
            (b'F', b'Y', 3), (b'F', b'W', 1),
            (b'P', b'A', -1),
            (b'S', b'T', 1),
            (b'Y', b'W', 2),
        ];
        
        for (aa1, aa2, score) in similar {
            m.insert((aa1, aa2), score);
            m.insert((aa2, aa1), score); // Symmetric
        }
        
        m
    };
}

/// Compute alignment score using BLOSUM62-like matrix
pub fn compute_alignment_score(aln: &Alignment) -> f64 {
    let query_bytes = aln.query.as_bytes();
    let target_bytes = aln.target.as_bytes();
    
    let mut score = 0i32;
    let mut qi = 0;
    let mut ti = 0;
    
    for op in &aln.operations {
        match op {
            EditOp::Match => {
                if qi < query_bytes.len() && ti < target_bytes.len() {
                    let aa1 = query_bytes[qi];
                    let aa2 = target_bytes[ti];
                    score += BLOSUM62.get(&(aa1, aa2)).copied().unwrap_or(-3);
                    qi += 1;
                    ti += 1;
                }
            }
            EditOp::Substitution => {
                if qi < query_bytes.len() && ti < target_bytes.len() {
                    let aa1 = query_bytes[qi];
                    let aa2 = target_bytes[ti];
                    score += BLOSUM62.get(&(aa1, aa2)).copied().unwrap_or(-3);
                    qi += 1;
                    ti += 1;
                }
            }
            EditOp::Insertion => {
                score -= 4; // Gap penalty
                ti += 1;
            }
            EditOp::Deletion => {
                score -= 4; // Gap penalty
                qi += 1;
            }
        }
    }
    
    score as f64
}

/// Compute normalized alignment score (0-1 range)
pub fn compute_normalized_score(aln: &Alignment) -> f64 {
    let raw_score = compute_alignment_score(aln);
    
    // Compute max possible score (perfect match)
    let query_bytes = aln.query.as_bytes();
    let mut max_score = 0i32;
    for &aa in query_bytes {
        max_score += BLOSUM62.get(&(aa, aa)).copied().unwrap_or(4);
    }
    
    let target_bytes = aln.target.as_bytes();
    let mut max_score_target = 0i32;
    for &aa in target_bytes {
        max_score_target += BLOSUM62.get(&(aa, aa)).copied().unwrap_or(4);
    }
    
    let max_possible = max_score.max(max_score_target);
    
    if max_possible == 0 {
        return 0.0;
    }
    
    // Normalize: negative scores -> 0, positive scores scale to 1
    let normalized = (raw_score - (aln.edit_distance as f64 * -4.0)) / (max_possible as f64);
    normalized.max(0.0).min(1.0)
}

/// Simple scoring: just count mismatches
pub fn simple_mismatch_score(aln: &Alignment) -> f64 {
    1.0 - (aln.edit_distance as f64 / aln.query.len().max(aln.target.len()) as f64)
}

/// Segment matching score
pub fn segment_match_score(query_segment: &str, db_segment: &str, normalize: bool) -> f64 {
    let query_norm = if normalize {
        query_segment.split('*').next().unwrap_or(query_segment)
    } else {
        query_segment
    };
    
    let db_norm = if normalize {
        db_segment.split('*').next().unwrap_or(db_segment)
    } else {
        db_segment
    };
    
    if query_norm == db_norm {
        1.0
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alignment::align;
    
    #[test]
    fn test_compute_alignment_score() {
        let aln = align("CASSLGQAYEQYF", "CASSLGQAYEQYF");
        let score = compute_alignment_score(&aln);
        assert!(score > 0.0);
        
        let aln = align("AAAA", "TTTT");
        let score = compute_alignment_score(&aln);
        assert!(score < 0.0);
    }
    
    #[test]
    fn test_segment_match_score() {
        assert_eq!(segment_match_score("TRBV12-3*01", "TRBV12-3*02", true), 1.0);
        assert_eq!(segment_match_score("TRBV12-3*01", "TRBV12-3*02", false), 0.0);
        assert_eq!(segment_match_score("TRBV12-3", "TRBV12-4", true), 0.0);
    }
}
