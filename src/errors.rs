//! All errors that can occur in the virolution library.

use std::fmt;

pub type Result<T> = std::result::Result<T, VirolutionError>;

#[derive(Clone, Debug)]
pub enum VirolutionError {
    ImplementationError(String),
    InitializationError(String),
    ReadError(String),
}

impl fmt::Display for VirolutionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            VirolutionError::ImplementationError(message) => {
                write!(f, "ImplementationError: {}", message)
            }
            VirolutionError::InitializationError(message) => {
                write!(f, "InitializationError: {}", message)
            }
            VirolutionError::ReadError(message) => {
                write!(f, "ReadError: {}", message)
            }
        }
    }
}

impl std::error::Error for VirolutionError {}
