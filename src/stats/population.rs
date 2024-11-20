use itertools::Itertools;

use crate::core::Population;
use crate::encoding::Symbol;

/// Trait extension to compute frequencies of nucleotides in a population
pub trait PopulationFrequencies<S: Symbol> {
    fn frequencies(&self) -> Vec<f64>;
}

impl<S: Symbol> PopulationFrequencies<S> for Population<S> {
    /// Compute the frequencies of each nucleotide in a population and return as a vector.
    fn frequencies(&self) -> Vec<f64> {
        let wildtype = &self[&0].get_wildtype().upgrade().unwrap();
        let sequence_length = wildtype.get_length();

        let population_size = self.len();

        let mut counts: Vec<i64> = vec![0; sequence_length * 4];

        for (haplotype_ref, haplotype_count) in self.iter().counts() {
            for (&pos, &change) in haplotype_ref.as_ref().get_mutations().iter() {
                counts[change as usize + S::SIZE * pos] += haplotype_count as i64;
            }
        }

        for (pos, symbol) in wildtype.get_sequence().iter().enumerate() {
            counts[symbol.index() + S::SIZE * pos] = population_size as i64
                - counts[S::SIZE * pos..S::SIZE * (pos + 1)]
                    .iter()
                    .sum::<i64>();
        }

        counts
            .iter()
            .map(|&x| x as f64 / population_size as f64)
            .collect()
    }
}

/// Trait extension to compute distances between to Populations
pub trait PopulationDistance<S: Symbol> {
    fn distance(&self, other: &Population<S>) -> f64;
}

impl<S: Symbol> PopulationDistance<S> for Population<S> {
    /// Compute the distance between this and another populations.
    ///
    /// The distances is the sum of absolutes of per-position nucleotide frequencies in each
    /// population.
    fn distance(&self, other: &Population<S>) -> f64 {
        self.frequencies()
            .iter()
            .zip(other.frequencies())
            .map(|(a, b)| (a - b).abs())
            .sum()
    }
}
