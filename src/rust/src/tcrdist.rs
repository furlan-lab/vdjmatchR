use serde::{Deserialize, Serialize};

/// BLOSUM62 substitution matrix for amino acid scoring
/// Amino acid order: A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
const BLOSUM62: [[i8; 20]; 20] = [
    // A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V
    [  4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0], // A
    [ -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3], // R
    [ -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3], // N
    [ -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3], // D
    [  0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1], // C
    [ -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2], // Q
    [ -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2], // E
    [  0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3], // G
    [ -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3], // H
    [ -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3], // I
    [ -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1], // L
    [ -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2], // K
    [ -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1], // M
    [ -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1], // F
    [ -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2], // P
    [  1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2], // S
    [  0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0], // T
    [ -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3], // W
    [ -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1], // Y
    [  0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4], // V
];

/// Convert amino acid character to BLOSUM62 matrix index
fn aa_to_index(aa: u8) -> Option<usize> {
    match aa {
        b'A' => Some(0),  b'R' => Some(1),  b'N' => Some(2),  b'D' => Some(3),
        b'C' => Some(4),  b'Q' => Some(5),  b'E' => Some(6),  b'G' => Some(7),
        b'H' => Some(8),  b'I' => Some(9),  b'L' => Some(10), b'K' => Some(11),
        b'M' => Some(12), b'F' => Some(13), b'P' => Some(14), b'S' => Some(15),
        b'T' => Some(16), b'W' => Some(17), b'Y' => Some(18), b'V' => Some(19),
        _ => None,
    }
}

/// Get BLOSUM62 score for two amino acids
fn blosum62_score(aa1: u8, aa2: u8) -> i8 {
    match (aa_to_index(aa1), aa_to_index(aa2)) {
        (Some(i1), Some(i2)) => BLOSUM62[i1][i2],
        _ => -4, // Default penalty for unknown amino acids
    }
}

/// Calculate tcrdist-style position score
/// Formula: max(0, 4 - BLOSUM62[a][b])
fn position_score(aa1: u8, aa2: u8) -> i32 {
    let blosum_score = blosum62_score(aa1, aa2) as i32;
    std::cmp::max(0, 4 - blosum_score)
}

/// Needleman-Wunsch alignment with tcrdist-style scoring
/// Returns alignment score (distance)
fn align_sequences(seq1: &str, seq2: &str, gap_penalty: i32) -> i32 {
    let seq1_bytes = seq1.as_bytes();
    let seq2_bytes = seq2.as_bytes();
    let len1 = seq1_bytes.len();
    let len2 = seq2_bytes.len();

    // Handle empty sequences
    if len1 == 0 && len2 == 0 {
        return 0;
    }
    if len1 == 0 {
        return (len2 as i32) * gap_penalty;
    }
    if len2 == 0 {
        return (len1 as i32) * gap_penalty;
    }

    // Initialize DP matrix
    let mut dp = vec![vec![0i32; len2 + 1]; len1 + 1];

    // Initialize first row and column
    for i in 0..=len1 {
        dp[i][0] = (i as i32) * gap_penalty;
    }
    for j in 0..=len2 {
        dp[0][j] = (j as i32) * gap_penalty;
    }

    // Fill DP matrix
    for i in 1..=len1 {
        for j in 1..=len2 {
            let match_score = dp[i - 1][j - 1] + position_score(seq1_bytes[i - 1], seq2_bytes[j - 1]);
            let delete_score = dp[i - 1][j] + gap_penalty;
            let insert_score = dp[i][j - 1] + gap_penalty;

            dp[i][j] = std::cmp::min(match_score, std::cmp::min(delete_score, insert_score));
        }
    }

    dp[len1][len2]
}

/// T-cell receptor with alpha and beta chain CDR sequences
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TCR {
    pub cdr1_a_aa: Option<String>,
    pub cdr2_a_aa: Option<String>,
    pub cdr3_a_aa: Option<String>,
    pub cdr1_b_aa: Option<String>,
    pub cdr2_b_aa: Option<String>,
    pub cdr3_b_aa: Option<String>,
}

impl TCR {
    pub fn new(
        cdr1_a_aa: Option<String>,
        cdr2_a_aa: Option<String>,
        cdr3_a_aa: Option<String>,
        cdr1_b_aa: Option<String>,
        cdr2_b_aa: Option<String>,
        cdr3_b_aa: Option<String>,
    ) -> Self {
        Self {
            cdr1_a_aa,
            cdr2_a_aa,
            cdr3_a_aa,
            cdr1_b_aa,
            cdr2_b_aa,
            cdr3_b_aa,
        }
    }
}

/// Calculate tcrdist distance between two TCRs
/// Combines alpha and beta chain distances
pub fn tcrdist(tcr1: &TCR, tcr2: &TCR) -> f64 {
    let alpha_dist = chain_distance(
        &tcr1.cdr1_a_aa,
        &tcr1.cdr2_a_aa,
        &tcr1.cdr3_a_aa,
        &tcr2.cdr1_a_aa,
        &tcr2.cdr2_a_aa,
        &tcr2.cdr3_a_aa,
    );

    let beta_dist = chain_distance(
        &tcr1.cdr1_b_aa,
        &tcr1.cdr2_b_aa,
        &tcr1.cdr3_b_aa,
        &tcr2.cdr1_b_aa,
        &tcr2.cdr2_b_aa,
        &tcr2.cdr3_b_aa,
    );

    (alpha_dist + beta_dist) as f64
}

/// Calculate distance for a single chain (alpha or beta)
fn chain_distance(
    cdr1_1: &Option<String>,
    cdr2_1: &Option<String>,
    cdr3_1: &Option<String>,
    cdr1_2: &Option<String>,
    cdr2_2: &Option<String>,
    cdr3_2: &Option<String>,
) -> i32 {
    let mut total_distance = 0;

    // CDR1 distance (weight = 1, gap penalty = 4)
    if let (Some(seq1), Some(seq2)) = (cdr1_1, cdr1_2) {
        total_distance += align_sequences(seq1, seq2, 4);
    }

    // CDR2 distance (weight = 1, gap penalty = 4)
    if let (Some(seq1), Some(seq2)) = (cdr2_1, cdr2_2) {
        total_distance += align_sequences(seq1, seq2, 4);
    }

    // CDR3 distance (weight = 3, gap penalty = 8)
    if let (Some(seq1), Some(seq2)) = (cdr3_1, cdr3_2) {
        total_distance += 3 * align_sequences(seq1, seq2, 8);
    }

    total_distance
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_blosum62_score() {
        // Test identical amino acids
        assert_eq!(blosum62_score(b'A', b'A'), 4);
        assert_eq!(blosum62_score(b'W', b'W'), 11);

        // Test different amino acids
        assert_eq!(blosum62_score(b'A', b'R'), -1);
        assert_eq!(blosum62_score(b'C', b'C'), 9);
    }

    #[test]
    fn test_position_score() {
        // Identical amino acids should have score 0
        assert_eq!(position_score(b'A', b'A'), 0);
        assert_eq!(position_score(b'W', b'W'), 0);

        // Different amino acids
        assert_eq!(position_score(b'A', b'R'), std::cmp::max(0, 4 - (-1)));
    }

    #[test]
    fn test_align_identical_sequences() {
        let score = align_sequences("CASSF", "CASSF", 4);
        assert_eq!(score, 0); // Identical sequences should have distance 0
    }

    #[test]
    fn test_align_different_sequences() {
        let score = align_sequences("CASS", "CASF", 4);
        assert!(score > 0); // Different sequences should have distance > 0
    }

    #[test]
    fn test_tcrdist_identical() {
        let tcr1 = TCR::new(
            Some("TGTGC".to_string()),
            Some("TGTGC".to_string()),
            Some("CASSF".to_string()),
            Some("TGTGC".to_string()),
            Some("TGTGC".to_string()),
            Some("CASSF".to_string()),
        );

        let dist = tcrdist(&tcr1, &tcr1);
        assert_eq!(dist, 0.0);
    }

    #[test]
    fn test_tcrdist_different() {
        let tcr1 = TCR::new(
            Some("TGTGC".to_string()),
            Some("TGTGC".to_string()),
            Some("CASSF".to_string()),
            Some("TGTGC".to_string()),
            Some("TGTGC".to_string()),
            Some("CASSF".to_string()),
        );

        let tcr2 = TCR::new(
            Some("TGTGA".to_string()),
            Some("TGTGA".to_string()),
            Some("CASSLF".to_string()),
            Some("TGTGA".to_string()),
            Some("TGTGA".to_string()),
            Some("CASSLF".to_string()),
        );

        let dist = tcrdist(&tcr1, &tcr2);
        assert!(dist > 0.0);
    }
}
