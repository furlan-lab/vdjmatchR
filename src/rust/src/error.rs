use thiserror::Error;

#[derive(Error, Debug)]
pub enum VdjMatchError {
    #[error("IO error: {0}")]
    Io(#[from] std::io::Error),
    
    #[error("CSV parsing error: {0}")]
    Csv(#[from] csv::Error),
    
    #[error("Database not found: {0}")]
    DatabaseNotFound(String),
    
    #[error("Invalid search scope: {0}")]
    InvalidSearchScope(String),
    
    #[error("Invalid filter expression: {0}")]
    InvalidFilter(String),
    
    #[error("Sequence error: {0}")]
    Sequence(String),
    
    #[error("Alignment error: {0}")]
    Alignment(String),
    
    #[error("Invalid configuration: {0}")]
    Configuration(String),
    
    #[error("Network error: {0}")]
    Network(#[from] reqwest::Error),
    
    #[error("Regex error: {0}")]
    Regex(#[from] regex::Error),
}

pub type Result<T> = std::result::Result<T, VdjMatchError>;
