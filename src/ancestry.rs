use std::cell::RefCell;
use std::collections::{HashMap, LinkedList};
use std::rc::Rc;
use std::rc::Weak;

struct Ancestry {
    ancestry: LinkedList<Ancestors>,
}

struct Ancestors {
    ancestors: Vec<usize>,
}

struct AncestorsTreeBuilderNode {
    children: Vec<RefCell<AncestorsTreeBuilderNode>>,
    pos: usize,
    dist: usize,
}

struct AncestorsTreeNode {
    parent: Option<Weak<AncestorsTreeNode>>,
    children: Vec<Rc<AncestorsTreeNode>>,
    id: usize,
    distance: usize,
}

trait DuplicateIndexMapExt {
    fn duplicate_index_map(&self) -> HashMap<usize, Vec<usize>>;
}

impl DuplicateIndexMapExt for Ancestors {
    fn duplicate_index_map(&self) -> HashMap<usize, Vec<usize>> {
        let mut indices: HashMap<usize, Vec<usize>> = HashMap::new();
        for (i, &v) in self.ancestors.iter().enumerate() {
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

    fn get_tree(&self) -> AncestorsTreeBuilderNode {
        // create a list of singular nodes
        let mut trees: LinkedList<RefCell<AncestorsTreeBuilderNode>> =
            (0..self.ancestry.back().unwrap().ancestors.len())
                .map(|pos| {
                    RefCell::new(AncestorsTreeBuilderNode {
                        children: Vec::new(),
                        pos,
                        dist: 0,
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
                            children: children,
                            pos,
                            dist: 1,
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
