//! Abstracted attribute computation and storage
//!
//! This module contains the abstractions necessary to compute and store attributes of a
//! `Haplotype`. The attributes are computed lazily and are thread safe. The attributes are
//! computed by `AttributeProvider` instances that are registered in an `AttributeSetDefinition`.
//! Attributes are stored in an `AttributeSet` instance that is derived from an
//! `AttributeSetDefinition` or other `AttributeSet` instances.

use std::borrow::Cow;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

use crate::encoding::Symbol;
use crate::errors::{Result, VirolutionError};
use crate::references::HaplotypeRef;

// Attribute value definition
#[derive(Clone, Debug)]
pub enum AttributeValue {
    U32(u32),
    U64(u64),
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
    String(String),
}

impl TryFrom<AttributeValue> for u32 {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<u32> {
        match value {
            AttributeValue::U32(value) => Ok(value),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to u32",
                value
            ))),
        }
    }
}

impl TryFrom<AttributeValue> for u64 {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<u64> {
        match value {
            AttributeValue::U64(value) => Ok(value),
            AttributeValue::U32(value) => Ok(value as u64),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to u32",
                value
            ))),
        }
    }
}

impl TryFrom<AttributeValue> for usize {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<usize> {
        match value {
            AttributeValue::U64(value) => Ok(value as usize),
            AttributeValue::U32(value) => Ok(value as usize),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to usize",
                value
            ))),
        }
    }
}

impl TryFrom<AttributeValue> for i32 {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<i32> {
        match value {
            AttributeValue::I32(value) => Ok(value),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to i32",
                value
            ))),
        }
    }
}

impl TryFrom<AttributeValue> for i64 {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<i64> {
        match value {
            AttributeValue::I64(value) => Ok(value),
            AttributeValue::I32(value) => Ok(value as i64),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to i64",
                value
            ))),
        }
    }
}

impl TryFrom<AttributeValue> for f32 {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<f32> {
        match value {
            AttributeValue::F32(value) => Ok(value),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to f32",
                value
            ))),
        }
    }
}

impl TryFrom<AttributeValue> for f64 {
    type Error = VirolutionError;

    fn try_from(value: AttributeValue) -> Result<f64> {
        match value {
            AttributeValue::F64(value) => Ok(value),
            AttributeValue::F32(value) => Ok(value as f64),
            _ => Err(VirolutionError::ValueError(format!(
                "Cannot convert attribute value ({:?}) to f64",
                value
            ))),
        }
    }
}

/// Definition of an attribute computation.
///
/// Any attribute provider must implement this trait to provide the attribute value in an attribute
/// set. Providers are collected in an attribute set definition (see `AttributeSetDefinition`) and
/// compute the attributes of a set (`AttributeSet`).
pub trait AttributeProvider<S: Symbol>: Sync + Send {
    fn compute(&self, haplotype: &HaplotypeRef<S>) -> AttributeValue;
}

/// Definition of an attribute set with providers.
///
/// The attribute set definition contains the providers that compute the attributes of the set.
/// It can be used to register new providers and to create an attribute set.
#[derive(Clone)]
pub struct AttributeSetDefinition<S: Symbol> {
    providers: HashMap<Cow<'static, str>, Arc<dyn AttributeProvider<S> + Send + Sync>>,
}

impl<S: Symbol> AttributeSetDefinition<S> {
    /// Create a new attribute set definition.
    pub fn new() -> Self {
        Self {
            providers: HashMap::new(),
        }
    }

    /// Register a new attribute provider.
    pub fn register(&mut self, name: &str, provider: Arc<dyn AttributeProvider<S> + Send + Sync>) {
        let identifier = Cow::Owned(name.to_string());
        self.providers.insert(identifier, provider);
    }

    /// Create a new attribute set.
    pub fn create(&self) -> AttributeSet<S> {
        AttributeSet {
            definition: Arc::new(self.clone()),
            values: RwLock::new(HashMap::new()),
        }
    }
}

/// Instance of an attribute set with storage.
///
/// Attribute sets can be derived from attribute set definitions and other attribute sets. They are
/// used to define and store the attributes of a `Haplotype`. The attributes are computed lazily
/// and are thread safe.
pub struct AttributeSet<S: Symbol> {
    definition: Arc<AttributeSetDefinition<S>>,
    values: RwLock<HashMap<Cow<'static, str>, AttributeValue>>,
}

impl<S: Symbol> AttributeSet<S> {
    /// Derive a new attribute set from an attribute set definition with new storage.
    pub fn derive(&self) -> AttributeSet<S> {
        AttributeSet {
            definition: self.definition.clone(),
            values: RwLock::new(HashMap::new()),
        }
    }

    /// Get or compute the value of an attribute.
    pub fn get_or_compute(&self, id: &str, haplotype: HaplotypeRef<S>) -> Result<AttributeValue> {
        let id_cow = Cow::Borrowed(id);

        // First, try to read the value
        {
            let values = self.values.read().unwrap();
            if let Some(value) = values.get(&id_cow) {
                return Ok(value.clone());
            }
        }

        // Compute the attribute using the provider
        if let Some(provider) = self.definition.providers.get(&id_cow) {
            let value = provider.compute(&haplotype);

            // Write the computed value
            let mut values = self.values.write().unwrap();
            values.insert(id_cow.into_owned().into(), value.clone());

            Ok(value)
        } else {
            // If there is no provider, return None
            Err(VirolutionError::ImplementationError(
                "No provider found for attribute".to_string(),
            ))
        }
    }

    /// Get the value of an attribute.
    pub fn get(&self, id: &str) -> Option<AttributeValue> {
        let id_cow = Cow::Borrowed(id);
        let values = self.values.read().unwrap();
        values.get(&id_cow).cloned()
    }
}

impl<S: Symbol> Clone for AttributeSet<S> {
    /// Clone an attribute set including the computed values.
    fn clone(&self) -> Self {
        let hash_map = self.values.read().unwrap();
        AttributeSet {
            definition: self.definition.clone(),
            values: RwLock::new(hash_map.clone()),
        }
    }
}

impl<S: Symbol> std::fmt::Debug for AttributeSet<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let values = self.values.read().unwrap();
        f.debug_struct("AttributeSet")
            .field("values", &values)
            .finish()
    }
}
