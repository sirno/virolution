use itertools::Itertools;

use crate::core::Population;
use crate::core::population::HaplotypeStore;
use crate::encoding::Symbol;
use crate::references::HaplotypeWeakTrait;

/// Trait extension to compute frequencies of nucleotides in a population
pub trait PopulationFrequencies {
    fn frequencies(&self) -> Vec<f64>;
}

impl<M: HaplotypeStore> PopulationFrequencies for Population<M> {
    /// Compute the frequencies of each nucleotide in a population and return as a vector.
    fn frequencies(&self) -> Vec<f64> {
        let wildtype = &self.get(&0).get_wildtype().upgrade().unwrap();
        let sequence_length = wildtype.get_length();

        let population_size = self.len();

        let mut counts: Vec<i64> = vec![0; sequence_length * 4];

        for (haplotype_ref, haplotype_count) in self.iter().counts() {
            for (&pos, &change) in haplotype_ref.get_mutations().iter() {
                counts[change as usize + M::Symbol::SIZE * pos] += haplotype_count as i64;
            }
        }

        for (pos, symbol) in wildtype.get_sequence().iter().enumerate() {
            counts[symbol.index() + M::Symbol::SIZE * pos] = population_size as i64
                - counts[M::Symbol::SIZE * pos..M::Symbol::SIZE * (pos + 1)]
                    .iter()
                    .sum::<i64>();
        }

        counts
            .iter()
            .map(|&x| x as f64 / population_size as f64)
            .collect()
    }
}

/// Trait extension to compute attribute statistics of a population
pub trait PopulationStatistics {
    fn mean(&self, name: &'static str) -> f64;
}

impl<M: HaplotypeStore> PopulationStatistics for Population<M> {
    /// Compute the mean of a given attribute in a population.
    fn mean(&self, name: &'static str) -> f64 {
        let mut sum = 0.0;
        for haplotype in self.iter() {
            sum += f64::try_from(haplotype.get_or_compute_attribute(name).unwrap()).unwrap();
        }
        sum / self.len() as f64
    }
}

/// Trait extension to compute distances between to Populations
pub trait PopulationDistance {
    fn distance(&self, other: &Self) -> f64;
}

impl<M: HaplotypeStore> PopulationDistance for Population<M> {
    /// Compute the distance between this and another populations.
    ///
    /// The distances is the sum of absolutes of per-position nucleotide frequencies in each
    /// population.
    fn distance(&self, other: &Self) -> f64 {
        self.frequencies()
            .iter()
            .zip(other.frequencies())
            .map(|(a, b)| (a - b).abs())
            .sum()
    }
}
