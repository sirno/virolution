//! Abstracted attribute computation and storage
//!
//! This module contains the abstractions necessary to compute and store attributes of a
//! `Haplotype`. The attributes are computed lazily and are thread safe. The attributes are
//! computed by `AttributeProvider` instances that are registered in an `AttributeSetDefinition`.
//! Attributes are stored in an `AttributeSet` instance that is derived from an
//! `AttributeSetDefinition` or other `AttributeSet` instances.

use derive_more::{Display, TryInto};
use std::collections::HashMap;
use std::path::Path;
use std::sync::{Arc, RwLock};

use crate::encoding::Symbol;
use crate::errors::{Result, VirolutionError};
use crate::references::{HaplotypeRef, HaplotypeWeak, HaplotypeWeakTrait};

// Attribute value definition
#[derive(Clone, Debug, TryInto, Display)]
pub enum AttributeValue {
    U32(u32),
    U64(u64),
    I32(i32),
    I64(i64),
    F32(f32),
    F64(f64),
    String(String),
}

/// Type of an attribute provider.
///
/// The type of an attribute provider determines when the attribute is computed. A lazy provider
/// computes the attribute only when it is requested. An eager provider computes the attribute
/// immediately when the attribute set is created.
pub enum AttributeProviderType {
    Lazy,
    Eager,
}

/// Definition of an attribute computation.
///
/// Any attribute provider must implement this trait to provide the attribute value in an attribute
/// set. Providers are collected in an attribute set definition (see `AttributeSetDefinition`) and
/// compute the attributes of a set (`AttributeSet`).
pub trait AttributeProvider<S: Symbol>: Sync + Send + std::fmt::Debug {
    /// Get the name of the attribute.
    ///
    /// This name should correspond to the key used to store the attribute in the attribute set.
    fn name(&self) -> &'static str;

    /// Get the type of the provider.
    ///
    /// The type determines when the attribute is computed. A lazy provider computes the attribute
    /// only when it is requested. An eager provider computes the attribute immediately when the
    /// attribute set is created.
    fn get_provider_type(&self) -> AttributeProviderType;

    /// Compute the attribute value for a haplotype.
    ///
    /// The haplotype is passed as an `Option` to allow for the computation of attributes that are
    /// not associated with a specific haplotype. In the case of a eager provider, the computation
    /// must not depend on the haplotype.
    fn compute(&self, haplotype: &Option<HaplotypeRef<S>>) -> AttributeValue;

    /// Optional mapping of the attribute value.
    fn map(&self, value: AttributeValue) -> AttributeValue {
        value
    }

    /// Store attribute provider in a file.
    // TODO: This method should be a separate trait.
    fn write(&self, path: &Path) -> Result<()>;
}

/// Definition of an attribute set with providers.
///
/// The attribute set definition contains the providers that compute the attributes of the set.
/// It can be used to register new providers and to create an attribute set.
#[derive(Clone, Debug)]
pub struct AttributeSetDefinition<S: Symbol> {
    providers: HashMap<&'static str, Arc<dyn AttributeProvider<S> + Send + Sync>>,
    eager: Vec<&'static str>,
}

impl<S: Symbol> Default for AttributeSetDefinition<S> {
    fn default() -> Self {
        Self::new()
    }
}

impl<S: Symbol> AttributeSetDefinition<S> {
    /// Create a new attribute set definition.
    pub fn new() -> Self {
        Self {
            providers: HashMap::new(),
            eager: Vec::new(),
        }
    }

    /// Get the providers of the definition.
    pub fn providers(&self) -> &HashMap<&'static str, Arc<dyn AttributeProvider<S> + Send + Sync>> {
        &self.providers
    }

    /// Register a new attribute provider.
    pub fn register(&mut self, provider: Arc<dyn AttributeProvider<S> + Send + Sync>) {
        // notify user if provider with the same name already exists and do not register it twice
        if self.providers.contains_key(provider.name()) {
            log::warn!("Provider with name {} already exists", provider.name());
            return;
        }

        match provider.get_provider_type() {
            AttributeProviderType::Lazy => self.register_lazy(provider),
            AttributeProviderType::Eager => self.register_eager(provider),
        }
    }

    /// Register a new eager attribute provider.
    fn register_eager(&mut self, provider: Arc<dyn AttributeProvider<S> + Send + Sync>) {
        let identifier = provider.name();
        self.providers.insert(identifier, provider);
        self.eager.push(identifier);
    }

    /// Register a new lazy attribute provider.
    fn register_lazy(&mut self, provider: Arc<dyn AttributeProvider<S> + Send + Sync>) {
        self.providers.insert(provider.name(), provider);
    }

    /// Create a new attribute set.
    pub fn create(&self, haplotype: HaplotypeWeak<S>) -> AttributeSet<S> {
        AttributeSet::new(Arc::new(self.clone()), haplotype)
    }
}

impl<S: Symbol> IntoIterator for AttributeSetDefinition<S> {
    type Item = Arc<dyn AttributeProvider<S> + Send + Sync>;
    type IntoIter = std::collections::hash_map::IntoValues<&'static str, Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.providers.into_values()
    }
}

impl<S: Symbol> Extend<Arc<dyn AttributeProvider<S> + Send + Sync>> for AttributeSetDefinition<S> {
    fn extend<T: IntoIterator<Item = Arc<dyn AttributeProvider<S> + Send + Sync>>>(
        &mut self,
        iter: T,
    ) {
        for provider in iter {
            self.register(provider);
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
    values: RwLock<HashMap<&'static str, AttributeValue>>,
    haplotype: HaplotypeWeak<S>,
}

impl<S: Symbol> AttributeSet<S> {
    fn new(
        definition: Arc<AttributeSetDefinition<S>>,
        haplotype: HaplotypeWeak<S>,
    ) -> AttributeSet<S> {
        let mut values = HashMap::new();

        for provider_id in definition.eager.iter() {
            let value = definition
                .providers
                .get(provider_id)
                .unwrap()
                .compute(&None);
            values.insert(*provider_id, value);
        }

        AttributeSet {
            definition,
            values: RwLock::new(values),
            haplotype,
        }
    }
    /// Derive a new attribute set from an attribute set definition with new storage.
    pub fn derive(&self, haplotype: HaplotypeWeak<S>) -> AttributeSet<S> {
        Self::new(self.definition.clone(), haplotype)
    }

    /// Get or compute the value of an attribute.
    pub fn get_or_compute_raw(&self, id: &'static str) -> Result<AttributeValue> {
        let provider = self.definition.providers.get(&id).ok_or_else(|| {
            VirolutionError::ImplementationError(format!("No provider found for attribute {id}"))
        })?;

        // First, try to read the value
        {
            let values = self.values.read().unwrap();
            if let Some(value) = values.get(&id) {
                return Ok(value.clone());
            }
        }

        // Compute the attribute using the provider
        let value = provider.compute(&self.haplotype.upgrade());

        // Write the computed value
        let mut values = self.values.write().unwrap();
        values.insert(id, value.clone());

        Ok(value)
    }

    /// Get or compute the value of an attribute.
    pub fn get_or_compute(&self, id: &'static str) -> Result<AttributeValue> {
        let provider = self.definition.providers.get(&id).ok_or_else(|| {
            VirolutionError::ImplementationError(format!("No provider found for attribute {id}"))
        })?;

        // First, try to read the value
        {
            let values = self.values.read().unwrap();
            if let Some(value) = values.get(&id) {
                return Ok(provider.map(value.clone()));
            }
        }

        // Compute the attribute using the provider
        let value = provider.compute(&self.haplotype.upgrade());

        // Write the computed value
        let mut values = self.values.write().unwrap();
        values.insert(id, value.clone());

        Ok(provider.map(value))
    }

    /// Get the value of an already computed attribute.
    ///
    /// Returns `None` if the attribute has not been computed yet, or if the attribute does not
    /// have an associated provider.
    pub fn get(&self, id: &str) -> Option<AttributeValue> {
        // Get the provider and return None if it doesn't exist.
        let provider = match self.definition.providers.get(id).ok_or_else(|| {
            VirolutionError::ImplementationError(format!("No provider found for attribute {id}"))
        }) {
            Ok(provider) => provider,
            Err(_) => return None,
        };
        let values = self.values.read().unwrap();
        values.get(&id).map(|value| provider.map(value.clone()))
    }

    /// Get a provider by name.
    pub fn get_provider(&self, id: &str) -> Option<Arc<dyn AttributeProvider<S> + Send + Sync>> {
        self.definition.providers.get(&id).cloned()
    }
}

impl<S: Symbol> Clone for AttributeSet<S> {
    /// Clone an attribute set including the computed values.
    fn clone(&self) -> Self {
        let hash_map = self.values.read().unwrap();
        AttributeSet {
            definition: self.definition.clone(),
            values: RwLock::new(hash_map.clone()),
            haplotype: self.haplotype.clone(),
        }
    }
}

impl<S: Symbol> std::fmt::Debug for AttributeSet<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let values = self.values.read().unwrap();
        let haplotype = self.haplotype.upgrade();
        f.debug_struct("AttributeSet")
            .field("haplotype", &haplotype)
            .field("values", &values)
            .finish()
    }
}
