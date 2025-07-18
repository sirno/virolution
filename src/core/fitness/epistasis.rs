#![allow(non_local_definitions, unexpected_cfgs)]
use npyz::WriterBuilder;
use smallvec::SmallVec;
use std::collections::HashMap;

use crate::core::haplotype::Changes;
use crate::encoding::Symbol;
use crate::errors::VirolutionError;
use crate::init::{FitnessDistribution, FitnessModel};
use crate::references::HaplotypeRef;

type EpistasisTableKey<S> = (usize, S);

#[derive(Debug, npyz::Deserialize, npyz::Serialize, npyz::AutoSerialize)]
pub struct EpiEntry {
    pos1: u64,
    base1: u8,
    pos2: u64,
    base2: u8,
    value: f64,
}

#[derive(Clone, Debug)]
pub struct EpistasisTable<S: Symbol> {
    table: HashMap<EpistasisTableKey<S>, HashMap<EpistasisTableKey<S>, f64>>,
}

impl<S: Symbol> EpistasisTable<S> {
    pub fn from_model(model: &FitnessModel) -> Result<Self, VirolutionError> {
        let table = match &model.distribution {
            FitnessDistribution::Epistatic(epi_params) => {
                let entries = epi_params.load_epistasis();
                Self::from_vec(entries)
            }
            _ => {
                return Err(VirolutionError::ImplementationError(
                    "Model is not epistatic".to_string(),
                ));
            }
        };
        Ok(table)
    }

    pub fn from_vec(entries: Vec<EpiEntry>) -> Self {
        let mut table: HashMap<EpistasisTableKey<S>, HashMap<EpistasisTableKey<S>, f64>> =
            HashMap::new();
        for entry in entries.iter() {
            table
                .entry((entry.pos1 as usize, S::decode(&entry.base1)))
                .or_default()
                .insert((entry.pos2 as usize, S::decode(&entry.base2)), entry.value);
            table
                .entry((entry.pos2 as usize, S::decode(&entry.base2)))
                .or_default()
                .insert((entry.pos1 as usize, S::decode(&entry.base1)), entry.value);
        }
        Self { table }
    }

    pub fn write(&self, writer: &mut dyn std::io::Write) -> Result<(), VirolutionError> {
        // Convert the table to a vector of entries
        let entries: Vec<EpiEntry> = self
            .table
            .iter()
            .flat_map(|((pos1, base1), interactions)| {
                interactions
                    .iter()
                    .filter_map(move |((pos2, base2), value)| {
                        if pos1 > pos2 {
                            None
                        } else {
                            Some(EpiEntry {
                                pos1: *pos1 as u64,
                                base1: base1.encode(),
                                pos2: *pos2 as u64,
                                base2: base2.encode(),
                                value: *value,
                            })
                        }
                    })
            })
            .collect();

        // Setup Writer
        let shape = vec![entries.len() as u64];
        let mut npy_writer = npyz::WriteOptions::new()
            .default_dtype()
            .shape(&shape)
            .writer(writer)
            .begin_nd()
            .map_err(|e| VirolutionError::WriteError(format!("{e}")))?;

        // Write entries to file
        npy_writer
            .extend(entries)
            .map_err(|e| VirolutionError::WriteError(format!("{e}")))?;
        npy_writer
            .finish()
            .map_err(|e| VirolutionError::WriteError(format!("{e}")))?;

        Ok(())
    }

    /// Compute factor to update the fitness of a haplotype change based on the epistasis table
    pub fn update_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        // Callers should ensure that the haplotype is mutant
        if !haplotype.is_mutant() {
            panic!("Cannot update fitness of haplotype that is not mutant");
        }

        // If there are not changes then there is no update
        let changes = match haplotype.try_get_changes() {
            Some(changes) => changes,
            None => {
                return 1.;
            }
        };

        // We evaluate this first, because if we do not need to evaluate any changes, we do not
        // have to call `get_mutations` on the haplotype which can be expensive.
        let candidate_changes = changes
            .iter()
            .filter(|changes| {
                self.table.contains_key(&(changes.position, changes.to))
                    || self.table.contains_key(&(changes.position, changes.to))
            })
            .cloned()
            .collect::<SmallVec<Changes<S>>>();

        // If no changes are present in the table, return the fitness
        if candidate_changes.is_empty() {
            return 1.;
        }

        // Collect all mutations of the haplotype
        let mutations = haplotype.get_mutations();
        let mut fitness = 1.;

        // Add any epistatic effects
        for change in &candidate_changes {
            // Get any interactions of the current mutation
            let interactions_add = self.table.get(&(change.position, change.to));
            let interactions_remove = self.table.get(&(change.position, change.from));

            // Apply any interactions
            for (pos, current) in mutations.iter() {
                // Deal with a rare case where two mutations have interactions with each
                // other. In this case, we enforce that the interaction is only applied
                // once, by enforcing an order when reading epistatic interactions.
                if pos <= &change.position && changes.iter().any(|c| c.position == *pos) {
                    continue;
                }

                if pos != &change.position {
                    if let Some(interaction) = interactions_add {
                        // If there is an interaction, multiply the fitness
                        if let Some(v) = interaction.get(&(*pos, S::decode(current))) {
                            fitness *= v;
                        }
                    }

                    if let Some(interaction) = interactions_remove {
                        // If there is an interaction, divide the fitness
                        if let Some(v) = interaction.get(&(*pos, S::decode(current))) {
                            fitness /= v;
                        }
                    }
                }
            }
        }

        fitness
    }

    /// Compute the fitness contribution within the epistasis table for a given haplotype
    pub fn compute_fitness(&self, haplotype: &HaplotypeRef<S>) -> f64 {
        let mut fitness = 1.;
        let mutations = haplotype.get_mutations();
        mutations.iter().for_each(|(position, current)| {
            if let Some(interaction) = self.table.get(&(*position, S::decode(current))) {
                interaction
                    .iter()
                    .filter(|((pos, base), _)| {
                        // Enforce that the interaction is only applied once
                        pos > position && mutations.get(pos) != Some(&(*base.index() as u8))
                    })
                    .for_each(|(_, value)| {
                        fitness *= value;
                    });
            }
        });
        fitness
    }
}

#[cfg(test)]
mod tests {}
