use crate::sequence::{Cdr3Sequence, SearchScope};
use std::cmp::min;

/// Edit distance and alignment operations
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EditOp {
    Match,
    Substitution,
    Insertion,
    Deletion,
}

#[derive(Debug, Clone)]
pub struct Alignment {
    pub query: String,
    pub target: String,
    pub operations: Vec<EditOp>,
    pub substitutions: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub edit_distance: usize,
}

impl Alignment {
    pub fn within_scope(&self, scope: &SearchScope) -> bool {
        self.substitutions <= scope.substitutions
            && self.insertions <= scope.insertions
            && self.deletions <= scope.deletions
            && self.edit_distance <= scope.total
    }
}

/// Compute edit distance between two sequences
pub fn edit_distance(seq1: &str, seq2: &str) -> usize {
    let len1 = seq1.len();
    let len2 = seq2.len();
    
    if len1 == 0 {
        return len2;
    }
    if len2 == 0 {
        return len1;
    }
    
    let seq1_bytes = seq1.as_bytes();
    let seq2_bytes = seq2.as_bytes();
    
    let mut prev_row: Vec<usize> = (0..=len2).collect();
    let mut curr_row: Vec<usize> = vec![0; len2 + 1];
    
    for i in 1..=len1 {
        curr_row[0] = i;
        
        for j in 1..=len2 {
            let cost = if seq1_bytes[i - 1] == seq2_bytes[j - 1] { 0 } else { 1 };
            
            curr_row[j] = min(
                min(
                    prev_row[j] + 1,      // deletion
                    curr_row[j - 1] + 1,  // insertion
                ),
                prev_row[j - 1] + cost,    // substitution/match
            );
        }
        
        std::mem::swap(&mut prev_row, &mut curr_row);
    }
    
    prev_row[len2]
}

/// Check if two sequences match within the given search scope using edit distance
pub fn matches_within_scope(query: &Cdr3Sequence, target: &Cdr3Sequence, scope: &SearchScope) -> bool {
    if scope.is_exact() {
        return query.sequence == target.sequence;
    }
    
    let distance = edit_distance(&query.sequence, &target.sequence);
    distance <= scope.total
}

/// Perform detailed alignment with operation tracking
pub fn align(query: &str, target: &str) -> Alignment {
    let len1 = query.len();
    let len2 = target.len();
    
    let query_bytes = query.as_bytes();
    let target_bytes = target.as_bytes();
    
    // Create DP table
    let mut dp = vec![vec![0usize; len2 + 1]; len1 + 1];
    
    // Initialize
    for i in 0..=len1 {
        dp[i][0] = i;
    }
    for j in 0..=len2 {
        dp[0][j] = j;
    }
    
    // Fill DP table
    for i in 1..=len1 {
        for j in 1..=len2 {
            let cost = if query_bytes[i - 1] == target_bytes[j - 1] { 0 } else { 1 };
            
            dp[i][j] = min(
                min(
                    dp[i - 1][j] + 1,      // deletion
                    dp[i][j - 1] + 1,      // insertion
                ),
                dp[i - 1][j - 1] + cost,   // substitution/match
            );
        }
    }
    
    // Backtrack to get operations
    let mut operations = Vec::new();
    let mut i = len1;
    let mut j = len2;
    let mut substitutions = 0;
    let mut insertions = 0;
    let mut deletions = 0;
    
    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let cost = if query_bytes[i - 1] == target_bytes[j - 1] { 0 } else { 1 };
            
            if dp[i][j] == dp[i - 1][j - 1] + cost {
                if cost == 0 {
                    operations.push(EditOp::Match);
                } else {
                    operations.push(EditOp::Substitution);
                    substitutions += 1;
                }
                i -= 1;
                j -= 1;
                continue;
            }
        }
        
        if i > 0 && dp[i][j] == dp[i - 1][j] + 1 {
            operations.push(EditOp::Deletion);
            deletions += 1;
            i -= 1;
        } else if j > 0 && dp[i][j] == dp[i][j - 1] + 1 {
            operations.push(EditOp::Insertion);
            insertions += 1;
            j -= 1;
        }
    }
    
    operations.reverse();
    
    Alignment {
        query: query.to_string(),
        target: target.to_string(),
        operations,
        substitutions,
        insertions,
        deletions,
        edit_distance: dp[len1][len2],
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_edit_distance() {
        assert_eq!(edit_distance("CASSLGQAYEQYF", "CASSLGQAYEQYF"), 0);
        assert_eq!(edit_distance("CASSLGQAYEQYF", "CASSLGQAYEQYY"), 1);
        assert_eq!(edit_distance("CASSLGQAYEQYF", "CASSGQAYEQYF"), 1);
        assert_eq!(edit_distance("", "ABC"), 3);
        assert_eq!(edit_distance("ABC", ""), 3);
    }
    
    #[test]
    fn test_matches_within_scope() {
        let seq1 = Cdr3Sequence::new("CASSLGQAYEQYF".to_string());
        let seq2 = Cdr3Sequence::new("CASSLGQAYEQYY".to_string());
        
        let scope = SearchScope { substitutions: 1, insertions: 0, deletions: 0, total: 1 };
        assert!(matches_within_scope(&seq1, &seq2, &scope));
        
        let scope = SearchScope::EXACT;
        assert!(!matches_within_scope(&seq1, &seq2, &scope));
    }
    
    #[test]
    fn test_align() {
        let aln = align("CASSLGQAYEQYF", "CASSLGQAYEQYY");
        assert_eq!(aln.substitutions, 1);
        assert_eq!(aln.insertions, 0);
        assert_eq!(aln.deletions, 0);
        assert_eq!(aln.edit_distance, 1);
    }
}
