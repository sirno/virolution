//! This module contains the Ancestry struct and its associated methods.
//!
//! The Ancestry struct is used to store the ancestors of a given node. It can
//! be used to reconstruct the tree of individual ancestors for a single
//! sampling time.

use std::cell::RefCell;
use std::collections::{HashMap, LinkedList};
use std::rc::Rc;
use std::rc::Weak;

type Ancestors = Vec<usize>;

pub struct Ancestry {
    ancestry: LinkedList<Ancestors>,
}

struct AncestorsTreeBuilderNode {
    pos: usize,
    dist: usize,
    name: Option<String>,
    children: Vec<RefCell<AncestorsTreeBuilderNode>>,
}

trait DuplicateIndexMapExt {
    fn duplicate_index_map(&self) -> HashMap<usize, Vec<usize>>;
}

impl DuplicateIndexMapExt for Ancestors {
    fn duplicate_index_map(&self) -> HashMap<usize, Vec<usize>> {
        let mut indices: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, &v) in self.iter().enumerate() {
            indices.entry(v).or_insert(vec![]).push(i);
        }
        indices
    }
}

impl Ancestry {
    fn new() -> Ancestry {
        Ancestry {
            ancestry: LinkedList::new(),
        }
    }

    fn add_ancestors(&mut self, ancestors: Ancestors) {
        self.ancestry.push_back(ancestors);
    }

    /// Creates a tree of ancestors from the ancestors with the leaf nodes given
    /// by the names.
    fn get_tree(&self, names: &[&str]) -> AncestorsTreeBuilderNode {
        // create a list of singular nodes
        let length = self.ancestry.back().unwrap().len();
        let mut trees: LinkedList<RefCell<AncestorsTreeBuilderNode>> = (0..length)
            .map(|pos| {
                RefCell::new(AncestorsTreeBuilderNode {
                    pos,
                    dist: 0,
                    name: Some(names[pos].to_string()),
                    children: Vec::new(),
                })
            })
            .collect();

        // continuously merge nodes until there is only one node left
        for ancestors in self.ancestry.iter().rev() {
            let index_map = ancestors.duplicate_index_map();

            for (pos, indices) in index_map {
                match indices {
                    // if there is only one index, then we can just update the node
                    // at that index
                    indices if indices.len() == 1 => {
                        let index = indices[0];
                        if let Some(node) = trees.iter().find(|node| node.borrow().pos == index) {
                            node.borrow_mut().dist += 1;
                            node.borrow_mut().pos = pos;
                        }
                    }
                    // if there is more than one index, then we need to merge the
                    // nodes at those indices
                    indices if indices.len() > 1 => {
                        let children: Vec<RefCell<AncestorsTreeBuilderNode>> = trees
                            .extract_if(|node| indices.contains(&node.borrow().pos))
                            .collect();
                        trees.push_back(RefCell::new(AncestorsTreeBuilderNode {
                            pos,
                            dist: 0,
                            name: None,
                            children: children,
                        }));
                    }
                    _ => unreachable!(),
                }
            }

            if trees.len() == 1 {
                break;
            }
        }
        trees.pop_back().unwrap().into_inner().into()
    }
}
