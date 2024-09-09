//! Caching for expensive operations

use cached::{Cached, TimedCache};
use derive_more::{Deref, DerefMut};
use std::sync::{Arc, Mutex};

#[derive(Debug, Deref, DerefMut)]
pub struct CachedValue<T>(pub Arc<T>);

impl<T> Clone for CachedValue<T> {
    fn clone(&self) -> Self {
        CachedValue(self.0.clone())
    }
}

impl<T> Default for CachedValue<T>
where
    T: Default,
{
    fn default() -> Self {
        CachedValue(Arc::new(Default::default()))
    }
}

impl<T> CachedValue<T> {
    pub fn new(value: T) -> Self {
        CachedValue(Arc::new(value))
    }
}

impl<T> CachedValue<T>
where
    T: std::fmt::Debug + std::clone::Clone,
{
    #[inline]
    pub fn clone_inner(self) -> T {
        Arc::unwrap_or_clone(self.0)
    }
}

pub struct VirolutionCache<T> {
    cache: Mutex<cached::TimedCache<usize, CachedValue<T>>>,
}

impl<T> VirolutionCache<T> {
    pub fn new(seconds: u64) -> Self {
        Self {
            cache: Mutex::new(TimedCache::with_lifespan_and_refresh(seconds, true)),
        }
    }

    pub fn cache_get(&self, key: &usize) -> Option<CachedValue<T>> {
        let mut cache = self.cache.lock().unwrap();
        cache.cache_get(key).cloned()
    }

    pub fn cache_set(&self, key: usize, value: T) -> CachedValue<T> {
        let cached_value = CachedValue::<T>(Arc::new(value));
        let mut cache = self.cache.lock().unwrap();
        cache.cache_set(key, cached_value.clone());
        cached_value
    }
}
