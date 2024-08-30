use crate::core::Population;
use itertools::Itertools;

/// Trait extension to compute frequencies of nucleotides in a population
pub trait PopulationFrequencies {
    fn frequencies(&self) -> Vec<f64>;
}

impl PopulationFrequencies for Population {
    /// Compute the frequencies of each nucleotide in a population and return as a vector.
    fn frequencies(&self) -> Vec<f64> {
        let wildtype = &self[&0].get_wildtype();
        let sequence_length = wildtype.get_length();

        let population_size = self.len();

        let mut counts: Vec<i64> = vec![0; sequence_length * 4];

        for (haplotype_ref, haplotype_count) in self.iter().counts() {
            for (&pos, &change) in haplotype_ref.as_ref().get_mutations().iter() {
                if let Some(symbol) = change {
                    counts[symbol as usize + 4 * pos] += haplotype_count as i64;
                }
            }
        }

        for (pos, symbol) in wildtype.get_sequence().iter().enumerate() {
            if let Some(symbol) = symbol {
                counts[*symbol as usize + 4 * pos] =
                    population_size as i64 - counts[4 * pos..4 * (pos + 1)].iter().sum::<i64>();
            }
        }

        counts
            .iter()
            .map(|&x| x as f64 / population_size as f64)
            .collect()
    }
}

/// Trait extension to compute distances between to Populations
pub trait PopulationDistance {
    fn distance(&self, other: &Population) -> f64;
}

impl PopulationDistance for Population {
    /// Compute the distance between this and another populations.
    ///
    /// The distances is the sum of absolutes of per-position nucleotide frequencies in each
    /// population.
    fn distance(&self, other: &Population) -> f64 {
        self.frequencies()
            .iter()
            .zip(other.frequencies())
            .map(|(a, b)| (a - b).abs())
            .sum()
    }
}
