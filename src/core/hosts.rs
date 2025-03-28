//! Host abstraction for the simulation.
//!
//! The host abstraction constitutes three core components:
//!
//! 1. `Host`: A trait that defines the host's interaction with infectants.
//! 2. `HostSpec`: A host specification that associates hosts with a specific `Host`
//!    implementation.
//! 3. `HostMap`: A map from host indices to infectant indices, which is used to efficiently
//!    invert the host-infectant relationship.
//!
use crate::encoding::Symbol;
use crate::references::HaplotypeRef;
use derive_more::{Deref, DerefMut};
use rand::prelude::*;
use std::ops::Range;

/// A host that can be infected by infectants.
///
/// This trait is used to define the individual host's effect on the infectants. Currently, this
/// is evaluated in three stages:
///
/// 1. Infection: Evaluate whether an infectant does infect the host.
/// 2. Mutation: Mutate multiple infectants that are currently infecting the host. This may include
///    recombinations.
/// 3. Replication: Compute the number of infectants that will be released from the host.
///
pub trait Host<S: Symbol>: std::fmt::Debug + Send + Sync + 'static {
    /// Decide whether the infectant should infect the host.
    fn infect(&self, haplotype: &HaplotypeRef<S>, rng: &mut ThreadRng) -> bool;

    /// Any mutations that the infectants should experience after
    fn mutate(&self, haplotype: &mut [HaplotypeRef<S>], rng: &mut ThreadRng);

    /// Compute the number of infectants that will be released.
    fn replicate(
        &self,
        haplotype: &[HaplotypeRef<S>],
        offspring: &mut [usize],
        rng: &mut ThreadRng,
    );

    /// Get the attributes that are associated with the host.
    fn get_attributes(&self) -> Vec<&'static str>;

    /// Clone the host.
    fn clone_box(&self) -> Box<dyn Host<S>>;
}

/// A collection of hosts that can be accessed by index.
#[derive(Clone, Debug, Deref, DerefMut)]
pub struct HostSpecs<S: Symbol>(Vec<HostSpec<S>>);

impl<S: Symbol> HostSpecs<S> {
    pub fn new() -> Self {
        Self(Vec::new())
    }

    pub fn from_vec(specs: Vec<HostSpec<S>>) -> Self {
        Self(specs)
    }

    pub fn try_get_spec_from_index(&self, index: usize) -> Option<&HostSpec<S>> {
        self.0.iter().find(|&spec| spec.range.contains(&index))
    }
}

impl<S: Symbol> Default for HostSpecs<S> {
    fn default() -> Self {
        Self::new()
    }
}

/// A host specification that associates a host with a specific `Host` implementation.
#[derive(Debug)]
pub struct HostSpec<S: Symbol> {
    pub range: Range<usize>,
    pub host: Box<dyn Host<S>>,
}

#[cfg(feature = "parallel")]
unsafe impl<S: Symbol> Send for HostSpec<S> {}
#[cfg(feature = "parallel")]
unsafe impl<S: Symbol> Sync for HostSpec<S> {}

impl<S: Symbol> HostSpec<S> {
    pub fn new(range: Range<usize>, host: Box<dyn Host<S>>) -> Self {
        Self { range, host }
    }
}

impl<S: Symbol> Clone for HostSpec<S> {
    fn clone(&self) -> Self {
        Self {
            range: self.range.clone(),
            host: self.host.clone_box(),
        }
    }
}

/// A buffer to store the host map in.
///
/// This buffer can be used to reuse memory for the host map between generations
pub struct HostMapBuffer {
    // number of hosts and infectants
    // this is used to set the limits of the buffer without resizing it
    n_hosts: usize,
    n_infectants: usize,

    // buffers
    infections: Vec<Option<usize>>,
    offsets: Vec<usize>,
    hosts: Vec<usize>,
}

impl std::fmt::Debug for HostMapBuffer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let n_infections = self.infections.iter().filter(|x| x.is_some()).count();
        let infectivity = n_infections as f64 / self.n_infectants as f64;
        f.debug_struct("HostMapBuffer")
            .field("n_hosts", &self.n_hosts)
            .field("n_infectants", &self.n_infectants)
            .field("n_infections", &n_infections)
            .field("infectivity", &infectivity)
            .finish()
    }
}

/// A map from host index to infectant indices.
///
/// The map is constructed from a list of infections between infectant and host indices. It uses
/// memory allocated by the caller to store the map, such that the map can reuse the memory, when
/// associations change.
impl HostMapBuffer {
    pub fn new(n_hosts: usize, n_infectants: usize) -> Self {
        Self {
            n_hosts,
            n_infectants,

            infections: vec![None; n_infectants],
            offsets: vec![0; n_hosts + 1],
            hosts: vec![0; n_infectants],
        }
    }

    /// Set the number of infectants that are currently active
    pub fn limit_infectants(&mut self, n_infectants: usize) {
        self.n_infectants = n_infectants;

        // resize the buffer if necessary
        if self.infections.len() < n_infectants {
            self.infections.resize(n_infectants, None);
            self.hosts.resize(n_infectants, 0);
        }
    }

    /// Set the number of hosts that are currently active
    pub fn limit_hosts(&mut self, n_hosts: usize) {
        self.n_hosts = n_hosts;

        // resize the buffer if necessary
        if self.offsets.len() < n_hosts + 1 {
            self.offsets.resize(n_hosts + 1, 0);
        }
    }

    pub fn build<F: FnMut((usize, &mut Option<usize>))>(&mut self, builder: F) {
        let infection_slice = &mut self.infections[0..self.n_infectants];

        // set new infections
        infection_slice.iter_mut().enumerate().for_each(builder);

        // reset offsets
        self.offsets[0..=self.n_hosts].fill(0);

        // compute prefix sum in place
        for host_idx in infection_slice.iter().flatten() {
            self.offsets[*host_idx] += 1;
        }

        for i in 0..self.n_hosts {
            self.offsets[i + 1] += self.offsets[i];
        }

        // fill host maps using backfill with offsets
        // after the loop, we will have offsets - count
        for (infectant_idx, maybe_host) in infection_slice.iter().enumerate() {
            if let Some(host_idx) = maybe_host {
                self.hosts[self.offsets[*host_idx] - 1] = infectant_idx;
                self.offsets[*host_idx] -= 1;
            }
        }
    }

    /// Get the slice of infectants that are associated with the host.
    pub fn get_slice(&self, host_idx: usize) -> &[usize] {
        if host_idx >= self.offsets.len() {
            return &[];
        }

        let start = self.offsets[host_idx];
        let end = self.offsets[host_idx + 1];

        &self.hosts[start..end]
    }

    /// Find the host that is associated with the infectant.
    pub fn find_host(&self, infectant_idx: usize) -> Option<usize> {
        self.infections[infectant_idx]
    }

    /// Iterate over infectants associated with each host
    pub fn iter_hosts(&self) -> impl Iterator<Item = &[usize]> + '_ {
        self.offsets
            .windows(2)
            .map(move |window| &self.hosts[window[0]..window[1]])
    }

    /// Iterate over the infectants associated with a range of hosts
    pub fn iter_range(&self, range: Range<usize>) -> impl Iterator<Item = &[usize]> + '_ {
        range.map(move |i| &self.hosts[self.offsets[i]..self.offsets[i + 1]])
    }
}

/// A map from host index to infectant indices.
///
/// The map is constructed from a list of infections between infectant and host indices. It uses
/// memory allocated by the caller to store the map, such that the map can reuse the memory, when
/// associations change.
pub struct HostMap<'a> {
    infections: &'a [Option<usize>],
    offsets: &'a [usize],
    hosts: &'a [usize],
}

impl<'a> HostMap<'a> {
    pub fn new(
        infections: &'a [Option<usize>],
        offsets: &'a mut [usize],
        hosts: &'a mut [usize],
    ) -> Self {
        // reset offsets
        offsets.fill(0);

        // compute prefix sum in place
        for host_idx in infections.iter().flatten() {
            offsets[*host_idx] += 1;
        }

        for i in 0..offsets.len() - 1 {
            offsets[i + 1] += offsets[i];
        }

        dbg!(&offsets);

        // fill host maps using backfill with offsets
        // after the loop, we will have offsets - count
        for (infectant_idx, maybe_host) in infections.iter().enumerate() {
            if let Some(host_idx) = maybe_host {
                dbg!(infectant_idx);
                hosts[offsets[*host_idx] - 1] = infectant_idx;
                offsets[*host_idx] -= 1;
            }
        }

        Self {
            infections,
            offsets,
            hosts,
        }
    }

    /// Get the slice of infectants that are associated with the host.
    pub fn get_slice(&self, host_idx: usize) -> &[usize] {
        if host_idx >= self.offsets.len() {
            return &[];
        }

        let start = self.offsets[host_idx];
        let end = self.offsets[host_idx + 1];

        &self.hosts[start..end]
    }

    /// Find the host that is associated with the infectant.
    pub fn find_host(&self, infectant_idx: usize) -> Option<usize> {
        self.infections[infectant_idx]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_host_map() {
        let associations: &[Option<usize>] = &[
            Some(0),
            Some(1),
            Some(0),
            None,
            Some(1),
            None,
            Some(0),
            Some(1),
            Some(0),
            Some(1),
        ];
        const N_HOSTS: usize = 2;
        const N_INFECTANTS: usize = 10;

        let mut offsets = [0; N_HOSTS + 1];
        let mut hosts = [0; N_INFECTANTS];

        let host_map = HostMap::new(associations, &mut offsets, &mut hosts);

        // test if we can get the correct slices
        assert_eq!(host_map.get_slice(0), &[8, 6, 2, 0]);
        assert_eq!(host_map.get_slice(1), &[9, 7, 4, 1]);

        // test if we can retrieve the correct hosts for each infectant
        assert_eq!(host_map.find_host(0), Some(0));
        assert_eq!(host_map.find_host(1), Some(1));
        assert_eq!(host_map.find_host(2), Some(0));
        assert_eq!(host_map.find_host(3), None);
        assert_eq!(host_map.find_host(4), Some(1));
        assert_eq!(host_map.find_host(5), None);
        assert_eq!(host_map.find_host(6), Some(0));
        assert_eq!(host_map.find_host(7), Some(1));
        assert_eq!(host_map.find_host(8), Some(0));
        assert_eq!(host_map.find_host(9), Some(1));
    }
}
