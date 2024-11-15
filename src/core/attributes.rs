//! AttributeSet module
//!
//! This module contains the `AttributeSet` struct which is used to store a set of attributes.
//! The `AttributeSet` is lazily initialized and is thread safe. It is used to store the attributes
//! of a `Haplotype` or a `Population` that are computed during the simulation.

use std::any::Any;
use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock, RwLock};

// Trait for attribute providers
trait AttributeProvider: Send + Sync {
    fn get_name(&self) -> &'static str;
    fn compute(&self, attributes: &AttributeSet) -> Arc<dyn Any + Send + Sync>;
}

// Static registry of attribute providers
fn get_attribute_providers(
) -> &'static OnceLock<RwLock<HashMap<&'static str, Arc<dyn AttributeProvider>>>> {
    static ATTRIBUTE_PROVIDERS: OnceLock<
        RwLock<HashMap<&'static str, Arc<dyn AttributeProvider>>>,
    > = OnceLock::new();
    &ATTRIBUTE_PROVIDERS
}

// AttributeSet struct
struct AttributeSet {
    attributes: Mutex<HashMap<&'static str, Arc<dyn Any + Send + Sync>>>,
}

impl AttributeSet {
    fn new() -> Self {
        AttributeSet {
            attributes: Mutex::new(HashMap::new()),
        }
    }

    fn get_attribute<T: 'static + Send + Sync>(&self, name: &'static str) -> Option<Arc<T>> {
        {
            // Check if attribute is already computed
            let attributes = self.attributes.lock().unwrap();
            if let Some(value) = attributes.get(name) {
                if let Ok(val) = value.clone().downcast::<T>() {
                    return Some(val);
                }
            }
        }

        // Attribute providers should be initialized when the static registry is accessed for the
        // first time
        assert!(get_attribute_providers().get().is_some());

        // Compute the attribute using the provider
        let providers = get_attribute_providers().get().unwrap().read().unwrap();
        if let Some(provider) = providers.get(name) {
            let value = provider.compute(self);
            // Store the computed value
            let mut attributes = self.attributes.lock().unwrap();
            attributes.insert(name, value.clone());

            if let Ok(val) = value.downcast::<T>() {
                return Some(val);
            }
        }

        None
    }
}

// Macro to register attribute providers
macro_rules! register_attribute_provider {
    ($provider:expr) => {{
        let mut providers = get_attribute_providers().get().unwrap().write().unwrap();
        providers.insert($provider.get_name(), Arc::new($provider));
    }};
}

// Example attribute provider
struct ExampleAttributeProvider;

impl AttributeProvider for ExampleAttributeProvider {
    fn get_name(&self) -> &'static str {
        "example_attribute"
    }

    fn compute(&self, _attributes: &AttributeSet) -> Arc<dyn Any + Send + Sync> {
        Arc::new(42u32)
    }
}

fn main() {
    // Register the attribute provider
    register_attribute_provider!(ExampleAttributeProvider);

    let attribute_set = AttributeSet::new();
    // Retrieve the attribute
    if let Some(value) = attribute_set.get_attribute::<u32>("example_attribute") {
        println!("Attribute value: {}", *value);
    } else {
        println!("Attribute not found");
    }
}
