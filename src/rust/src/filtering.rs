use crate::database::DatabaseEntry;
use regex::Regex;

/// Text filter for database columns
pub trait TextFilter {
    fn matches(&self, entry: &DatabaseEntry) -> bool;
}

/// Exact text match filter
pub struct ExactFilter {
    pub column: String,
    pub value: String,
}

impl TextFilter for ExactFilter {
    fn matches(&self, entry: &DatabaseEntry) -> bool {
        match self.column.as_str() {
            "species" => entry.species.eq_ignore_ascii_case(&self.value),
            "gene" => entry.gene.eq_ignore_ascii_case(&self.value),
            "antigen.species" => entry.antigen_species.eq_ignore_ascii_case(&self.value),
            "antigen.epitope" => entry.antigen_epitope.eq_ignore_ascii_case(&self.value),
            _ => false,
        }
    }
}

/// Regex filter for database columns
pub struct RegexFilter {
    pub column: String,
    pub pattern: Regex,
}

impl TextFilter for RegexFilter {
    fn matches(&self, entry: &DatabaseEntry) -> bool {
        let text = match self.column.as_str() {
            "species" => &entry.species,
            "gene" => &entry.gene,
            "antigen.species" => &entry.antigen_species,
            "antigen.epitope" => &entry.antigen_epitope,
            "antigen.gene" => {
                if let Some(ref gene) = entry.antigen_gene {
                    gene
                } else {
                    return false;
                }
            }
            _ => return false,
        };
        
        self.pattern.is_match(text)
    }
}

/// Parse and apply filter expression
/// Format: "__column__=~'pattern'" or "__column__=='value'"
pub fn parse_filter_expression(expr: &str) -> Result<Box<dyn TextFilter>, String> {
    // Simple parser for filter expressions
    // Supports: __column__=~'regex' and __column__=='value'
    
    if let Some(regex_match) = Regex::new(r"__([^_]+)__=~'([^']+)'").ok() {
        if let Some(captures) = regex_match.captures(expr) {
            let column = captures.get(1).unwrap().as_str().to_string();
            let pattern_str = captures.get(2).unwrap().as_str();
            
            if let Ok(pattern) = Regex::new(pattern_str) {
                return Ok(Box::new(RegexFilter { column, pattern }));
            }
        }
    }
    
    if let Some(exact_match) = Regex::new(r"__([^_]+)__=='([^']+)'").ok() {
        if let Some(captures) = exact_match.captures(expr) {
            let column = captures.get(1).unwrap().as_str().to_string();
            let value = captures.get(2).unwrap().as_str().to_string();
            
            return Ok(Box::new(ExactFilter { column, value }));
        }
    }
    
    Err(format!("Invalid filter expression: {}", expr))
}
