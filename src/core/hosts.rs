//! Host abstraction for the simulation.
//!
//! The host abstraction constitutes three core components:
//!
//! 1. `Host`: A trait that defines the host's interaction with infectants.
//! 2. `HostSpec`: A host specification that associates hosts with a specific `Host`
//!    implementation.
//! 3. `HostMap`: A map from host indices to infectant indices, which is used to efficiently
//!    invert the host-infectant relationship.
use crate::encoding::{Nucleotide, Symbol};
use crate::references::HaplotypeRef;

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
pub trait Host<S: Symbol> {
    /// Decide whether the infectant should infect the host.
    fn infect(&self, haplotype: &HaplotypeRef<S>) -> bool;

    /// Any mutations that the infectants should experience after
    fn mutate(&self, haplotype: &mut [HaplotypeRef<S>]);

    /// Compute the number of infectants that will be released.
    fn replicate(&self, haplotype: &[HaplotypeRef<S>], offspring: &mut [usize]);
}

/// A host specification that associates a host with a specific `Host` implementation.
pub struct HostSpec<S: Symbol = Nucleotide> {
    range: Range<usize>,
    host: dyn Host<S>,
}

/// A map from host index to infectant indices.
///
/// The map is constructed from a list of associations between infectant and host indices. It uses
/// memory allocated by the caller to store the map, such that the map can reuse the memory, when
/// associations change.
pub struct HostMap<'a> {
    associations: &'a [Option<usize>],
    offsets: &'a [usize],
    hosts: &'a [usize],
}

impl<'a> HostMap<'a> {
    pub fn new(
        associations: &'a [Option<usize>],
        offsets: &'a mut [usize],
        hosts: &'a mut [usize],
    ) -> Self {
        // reset offsets
        offsets.fill(0);

        // compute prefix sum in place
        for host_idx in associations.iter().flatten() {
            offsets[*host_idx] += 1;
        }

        for i in 0..offsets.len() - 1 {
            offsets[i + 1] += offsets[i];
        }

        dbg!(&offsets);

        // fill host maps using backfill with offsets
        // after the loop, we will have offsets - count
        for (infectant_idx, maybe_host) in associations.iter().enumerate() {
            if let Some(host_idx) = maybe_host {
                dbg!(infectant_idx);
                hosts[offsets[*host_idx] - 1] = infectant_idx;
                offsets[*host_idx] -= 1;
            }
        }

        Self {
            associations,
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
        self.associations[infectant_idx]
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
