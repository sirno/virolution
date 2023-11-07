//! This module contains the Ancestry struct and its associated methods.
//!
//! The Ancestry struct is used to store the ancestors of a given node. It can
//! be used to reconstruct the tree of individual ancestors for a single
//! sampling time.

use std::cell::RefCell;
use std::collections::{HashMap, LinkedList};

type Ancestors = Vec<usize>;

#[derive(Debug)]
pub struct Ancestry {
    ancestry: LinkedList<Ancestors>,
}

#[derive(Clone, Debug)]
struct AncestorsTreeNode {
    pos: usize,
    dist: usize,
    name: Option<String>,
    children: Vec<RefCell<AncestorsTreeNode>>,
}

trait DuplicateIndexMapExt {
    fn duplicate_index_map(&self) -> HashMap<usize, Vec<usize>>;
}

impl DuplicateIndexMapExt for Ancestors {
    fn duplicate_index_map(&self) -> HashMap<usize, Vec<usize>> {
        let mut indices: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, &v) in self.iter().enumerate() {
            indices.entry(v).or_default().push(i);
        }
        indices
    }
}

impl Default for Ancestry {
    fn default() -> Self {
        Self::new()
    }
}

impl Ancestry {
    pub fn new() -> Ancestry {
        Ancestry {
            ancestry: LinkedList::new(),
        }
    }

    /// Adds a new set of ancestors to the ancestry.
    ///
    /// Each set of ancestors is a vector of indices of the ancestors of a
    /// previous generation.
    pub fn add_ancestors(&mut self, ancestors: Ancestors) {
        self.ancestry.push_back(ancestors);
    }

    /// Creates a tree of ancestors from the ancestors with the leaf nodes given
    /// by the names.
    pub fn get_tree<T: AsRef<str>>(&self, names: &[T]) -> String {
        // create a list of singular nodes
        let length = self.ancestry.back().unwrap().len();
        let mut trees: LinkedList<RefCell<AncestorsTreeNode>> = (0..length)
            .map(|pos| {
                RefCell::new(AncestorsTreeNode {
                    pos,
                    dist: 0,
                    name: Some(names[pos].as_ref().to_string()),
                    children: Vec::new(),
                })
            })
            .collect();

        // continuously merge nodes until there is only one node left
        for ancestors in self.ancestry.iter().rev() {
            let index_map = ancestors.duplicate_index_map();
            let mut next: LinkedList<RefCell<AncestorsTreeNode>> = LinkedList::new();

            trees.iter().for_each(|node| {
                node.borrow_mut().dist += 1;
            });

            for (pos, indices) in index_map {
                match indices {
                    // if there is only one index, then we can just update the node
                    // at that index
                    indices if indices.len() == 1 => {
                        let index = indices[0];
                        if let Some(node) = trees
                            .extract_if(|node| node.borrow().pos == index)
                            .collect::<LinkedList<_>>()
                            .pop_back()
                        {
                            node.borrow_mut().pos = pos;
                            next.push_back(node);
                        }
                    }
                    // if there is more than one index, then we need to merge the
                    // nodes at those indices
                    indices if indices.len() > 1 => {
                        let mut children: LinkedList<RefCell<AncestorsTreeNode>> = trees
                            .extract_if(|node| indices.contains(&node.borrow().pos))
                            .collect();
                        match children.len() {
                            0 => {}
                            1 => {
                                let child = children.pop_back().unwrap();
                                child.borrow_mut().pos = pos;
                                next.push_back(child);
                            }
                            _ => {
                                let node = RefCell::new(AncestorsTreeNode {
                                    pos,
                                    dist: 0,
                                    name: None,
                                    children: children.into_iter().collect(),
                                });
                                next.push_back(node);
                            }
                        };
                    }
                    _ => unreachable!(),
                }
            }

            trees = next;

            if trees.len() == 1 {
                break;
            }
        }

        // merge remaining trees
        let node = match trees.len() {
            0 => panic!("No trees found"),
            1 => trees.pop_back().unwrap(),
            _ => {
                let children: Vec<RefCell<AncestorsTreeNode>> = trees.into_iter().collect();
                RefCell::new(AncestorsTreeNode {
                    pos: 0,
                    dist: 0,
                    name: None,
                    children,
                })
            }
        };

        format!("{};", node.borrow().to_newick())
    }
}

impl std::fmt::Display for AncestorsTreeNode {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(fmt, "{};", self.to_newick())
    }
}

impl AncestorsTreeNode {
    fn to_newick(&self) -> String {
        let mut newick = String::new();
        if !self.children.is_empty() {
            newick.push('(');
            for (i, child) in self.children.iter().enumerate() {
                if i > 0 {
                    newick.push(',');
                }
                newick.push_str(&child.borrow().to_newick());
            }
            newick.push(')');
        }
        if self.name.is_some() {
            newick.push_str(self.name.as_ref().unwrap());
        }
        if self.dist > 0 {
            newick.push_str(&format!(":{}", self.dist));
        }
        newick
    }
}
