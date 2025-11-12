use serde::{Deserialize, Serialize};
use std::fmt;

/// Represents a CDR3 amino acid sequence
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Cdr3Sequence {
    pub sequence: String,
}

impl Cdr3Sequence {
    pub fn new(sequence: String) -> Self {
        Self { sequence: sequence.to_uppercase() }
    }
    
    pub fn len(&self) -> usize {
        self.sequence.len()
    }
    
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
    
    pub fn as_bytes(&self) -> &[u8] {
        self.sequence.as_bytes()
    }
}

impl fmt::Display for Cdr3Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.sequence)
    }
}

/// Represents a T-cell receptor clonotype
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Clonotype {
    pub count: usize,
    pub frequency: f64,
    pub cdr3_aa: Cdr3Sequence,
    pub cdr3_nt: Option<String>,
    pub v_segment: String,
    pub d_segment: Option<String>,
    pub j_segment: String,
    pub sample_id: Option<String>,
    pub id_in_sample: Option<usize>,
}

impl Clonotype {
    pub fn new(
        cdr3_aa: String,
        v_segment: String,
        j_segment: String,
        count: usize,
        frequency: f64,
    ) -> Self {
        Self {
            count,
            frequency,
            cdr3_aa: Cdr3Sequence::new(cdr3_aa),
            cdr3_nt: None,
            v_segment,
            d_segment: None,
            j_segment,
            sample_id: None,
            id_in_sample: None,
        }
    }
    
    /// Normalize segment names (remove allele information)
    pub fn normalize_segment(segment: &str) -> String {
        segment.split('*').next().unwrap_or(segment).to_string()
    }
    
    pub fn v_normalized(&self) -> String {
        Self::normalize_segment(&self.v_segment)
    }
    
    pub fn j_normalized(&self) -> String {
        Self::normalize_segment(&self.j_segment)
    }
}

/// Represents search scope parameters for fuzzy matching
#[derive(Debug, Clone, Copy)]
pub struct SearchScope {
    pub substitutions: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub total: usize,
}

impl SearchScope {
    pub const EXACT: Self = Self {
        substitutions: 0,
        insertions: 0,
        deletions: 0,
        total: 0,
    };
    
    /// Parse search scope from string like "2,1,2,3" (s,i,d,t) or "2,2,3" (s,id,t)
    pub fn parse(s: &str) -> Result<Self, String> {
        let parts: Vec<&str> = s.split(',').collect();
        
        match parts.len() {
            3 => {
                // Format: s,id,t
                let substitutions = parts[0].parse().map_err(|_| format!("Invalid substitutions: {}", parts[0]))?;
                let indels = parts[1].parse().map_err(|_| format!("Invalid indels: {}", parts[1]))?;
                let total = parts[2].parse().map_err(|_| format!("Invalid total: {}", parts[2]))?;
                
                Ok(Self {
                    substitutions,
                    insertions: indels,
                    deletions: indels,
                    total,
                })
            }
            4 => {
                // Format: s,i,d,t
                let substitutions = parts[0].parse().map_err(|_| format!("Invalid substitutions: {}", parts[0]))?;
                let insertions = parts[1].parse().map_err(|_| format!("Invalid insertions: {}", parts[1]))?;
                let deletions = parts[2].parse().map_err(|_| format!("Invalid deletions: {}", parts[2]))?;
                let total = parts[3].parse().map_err(|_| format!("Invalid total: {}", parts[3]))?;
                
                Ok(Self {
                    substitutions,
                    insertions,
                    deletions,
                    total,
                })
            }
            _ => Err(format!("Invalid search scope format: {}", s)),
        }
    }
    
    pub fn is_exact(&self) -> bool {
        self.total == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_search_scope_parse() {
        let scope = SearchScope::parse("2,1,2,3").unwrap();
        assert_eq!(scope.substitutions, 2);
        assert_eq!(scope.insertions, 1);
        assert_eq!(scope.deletions, 2);
        assert_eq!(scope.total, 3);
        
        let scope = SearchScope::parse("2,2,3").unwrap();
        assert_eq!(scope.substitutions, 2);
        assert_eq!(scope.insertions, 2);
        assert_eq!(scope.deletions, 2);
        assert_eq!(scope.total, 3);
    }
}
