//! Error types and implementations

use std::fmt;

pub type Result<T> = std::result::Result<T, VirolutionError>;

#[derive(Clone, Debug)]
pub enum VirolutionError {
    ImplementationError(String),
    InitializationError(String),
    ValueError(String),
    VariantMissmatch(String),
    ReadError(String),
    WriteError(String),
}

impl fmt::Display for VirolutionError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            VirolutionError::ImplementationError(message) => {
                write!(f, "ImplementationError: {message}")
            }
            VirolutionError::InitializationError(message) => {
                write!(f, "InitializationError: {message}")
            }
            VirolutionError::ValueError(message) => {
                write!(f, "ValueError: {message}")
            }
            VirolutionError::VariantMissmatch(message) => {
                write!(f, "VariantMissmatch: {message}")
            }
            VirolutionError::ReadError(message) => {
                write!(f, "ReadError: {message}")
            }
            VirolutionError::WriteError(message) => {
                write!(f, "WriteError: {message}")
            }
        }
    }
}

impl std::error::Error for VirolutionError {}
