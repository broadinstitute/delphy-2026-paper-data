use itertools::EitherOrBoth;
use rand::Rng;
use serde::Serialize;
use std::collections::HashMap;

use crate::newick::NewickTree;
use crate::trees::{NodeLike, TraversalAction, TreeLike};

/// A clade fingerprint represents a set of tips.
/// Singletons are represented directly as a random fingerprint.
/// Larger sets are represented as the XOR of all of their tips' fingerprints.
#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Copy, Serialize)]
pub struct CladeFp(u64);

impl CladeFp {
    pub fn new(value: u64) -> CladeFp {
        CladeFp(value)
    }

    pub fn empty() -> CladeFp {
        CladeFp(0)
    }

    pub fn random(rng: &mut dyn Rng) -> CladeFp {
        CladeFp(rng.next_u64())
    }

    pub fn union(&self, other: &CladeFp) -> CladeFp {
        CladeFp(self.0 ^ other.0)
    }

    pub fn is_empty(&self) -> bool {
        self.0 == 0
    }
}

pub fn assign_tip_fps(tree: &NewickTree, rng: &mut dyn Rng) -> HashMap<String, CladeFp> {
    let mut result = HashMap::new();
    for node_ref in tree.any_order_iter() {
        let node = node_ref.borrow();
        if node.is_tip() {
            assert!(!node.name.is_empty(), "Tip with no name!");
            result.insert(node.name.clone(), CladeFp::random(rng));
        }
    }
    result
}

/// Compute the Robinson-Foulds distance between two trees represented as
/// sorted slices of clade fingerprints.
pub fn calc_rf_dist(sorted_split_fps_a: &[CladeFp], sorted_split_fps_b: &[CladeFp]) -> u64 {
    use EitherOrBoth::{Both, Left, Right};
    let mut distance = 0;
    for pair in itertools::merge_join_by(
        sorted_split_fps_a,
        sorted_split_fps_b,
        |a, b| a.cmp(b),
    ) {
        match pair {
            Both(_, _) => (),
            Left(_) | Right(_) => distance += 1,
        }
    }
    distance
}

#[derive(Debug, Serialize, Clone, Eq, PartialEq)]
pub enum CladeDefinition {
    TipClade {
        name: String,
        fp: CladeFp,
    },
    InnerNodeClade {
        subclade1: CladeFp,
        subclade2: CladeFp,
    },
}

impl CladeDefinition {
    pub fn fp(&self) -> CladeFp {
        match self {
            CladeDefinition::TipClade { fp, .. } => *fp,
            CladeDefinition::InnerNodeClade {
                subclade1,
                subclade2,
                ..
            } => subclade1.union(subclade2),
        }
    }

    fn is_complete(&self) -> bool {
        match self {
            CladeDefinition::TipClade { .. } => true,
            CladeDefinition::InnerNodeClade {
                subclade1,
                subclade2,
                ..
            } => !subclade1.is_empty() && !subclade2.is_empty(),
        }
    }
}

pub type CladeMap = HashMap<CladeFp, CladeDefinition>;

pub fn tips_in_clade(fp: &CladeFp, clade_map: &CladeMap) -> Vec<String> {

    fn go(fp: &CladeFp, clade_map: &CladeMap, result: &mut Vec<String>) {
        match clade_map.get(fp) {
            None => { },
            Some(CladeDefinition::TipClade { name, .. }) => result.push(name.clone()),
            Some(CladeDefinition::InnerNodeClade {
                subclade1,
                subclade2,
                ..
            }) => {
                go(subclade1, clade_map, result);
                go(subclade2, clade_map, result);
            }
        }
    }

    let mut result: Vec<String> = Vec::new();
    go(fp, clade_map, &mut result);
    result.sort();
    result
}

pub struct CladeInTreeInfo {
    pub time_from_root: f64,
}

/// Analyze the clades in a tree, populating a shared `CladeMap` with clade
/// definitions (only adding clades not already present) and returning per-tree
/// information (time from root) for each clade.
pub fn analyze_tree_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    clade_map: &mut CladeMap,
) -> HashMap<CladeFp, CladeInTreeInfo> {

    let mut node_clade_defn = CladeDefinition::InnerNodeClade {
        subclade1: CladeFp::empty(),
        subclade2: CladeFp::empty(),
    };
    let mut node_clade_defn_stack: Vec<CladeDefinition> = vec![];

    let mut time_from_root = 0.0;
    let mut time_from_root_stack: Vec<f64> = vec![];

    let mut per_tree_info: HashMap<CladeFp, CladeInTreeInfo> = HashMap::new();

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                // Shift focus from parent to node

                time_from_root_stack.push(time_from_root);
                if node_ref != tree.root {
                    // root-to-root distance is 0 by definition; ignore branch length then
                    time_from_root += node.branch_length
                };

                let new_clade_defn = {
                    if node.is_tip() {
                        CladeDefinition::TipClade {
                            name: node.name.clone(),
                            fp: *tip_fps
                                .get(node.name.as_str())
                                .expect("Unknown tip name"),
                        }
                    } else {
                        CladeDefinition::InnerNodeClade {
                            subclade1: CladeFp::empty(),
                            subclade2: CladeFp::empty(),
                        }
                    }
                };
                node_clade_defn_stack.push(std::mem::replace(&mut node_clade_defn, new_clade_defn));

                // Invariant: at this point, time_from_root and node_clade_defn refer to the
                // time and (partial) fingerprint for the current node after entering it
            }
            TraversalAction::Exit => {
                // Invariant: at this point, time_from_root and node_clade_defn refer to the
                // time and (partial) fingerprint for the current node before exiting it

                assert!(node_clade_defn.is_complete());
                let node_fp = node_clade_defn.fp();

                per_tree_info.insert(node_fp, CladeInTreeInfo { time_from_root });

                // After this line, node_clade_defn becomes the parent's clade definition
                let child_clade_defn =
                    std::mem::replace(&mut node_clade_defn, node_clade_defn_stack.pop().unwrap());
                clade_map.entry(node_fp).or_insert(child_clade_defn);

                // Shift focus from node to parent, merging this node's fingerprint into it

                // node_clade_defn is already the parent's clade definition; merge exited child's
                // fingerprint if necessary
                node_clade_defn = match node_clade_defn {
                    CladeDefinition::TipClade { .. } => node_clade_defn,
                    CladeDefinition::InnerNodeClade {
                        subclade1,
                        subclade2,
                    } => {
                        if subclade1.is_empty() {
                            CladeDefinition::InnerNodeClade {
                                subclade1: node_fp,
                                subclade2: CladeFp::empty(),
                            }
                        } else if subclade2.is_empty() {
                            CladeDefinition::InnerNodeClade {
                                subclade1,
                                subclade2: node_fp,
                            }
                        } else {
                            panic!("Non-bifurcating tree?")
                        }
                    }
                };

                time_from_root = time_from_root_stack.pop().unwrap();
            }
        }
    }

    per_tree_info
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand_pcg::Pcg64Mcg;
    use std::collections::HashSet;

    #[test]
    fn random_produces_distinct_nonempty_values() {
        let mut rng = Pcg64Mcg::seed_from_u64(42);
        let fps: Vec<CladeFp> = (0..100).map(|_| CladeFp::random(&mut rng)).collect();
        for fp in &fps {
            assert!(!fp.is_empty());
        }
        let unique: HashSet<CladeFp> = fps.into_iter().collect();
        assert_eq!(unique.len(), 100);
    }

    #[test]
    fn union_is_commutative() {
        let mut rng = Pcg64Mcg::seed_from_u64(1);
        let a = CladeFp::random(&mut rng);
        let b = CladeFp::random(&mut rng);
        assert_eq!(a.union(&b), b.union(&a));
    }

    #[test]
    fn union_is_associative() {
        let mut rng = Pcg64Mcg::seed_from_u64(1);
        let a = CladeFp::random(&mut rng);
        let b = CladeFp::random(&mut rng);
        let c = CladeFp::random(&mut rng);
        assert_eq!(a.union(&b).union(&c), a.union(&b.union(&c)));
    }

    #[test]
    fn subset_fingerprints_are_unique() {
        // Given N random tip fingerprints, every non-empty subset should produce
        // a distinct union fingerprint (with overwhelming probability).
        let mut rng = Pcg64Mcg::seed_from_u64(123);
        let n = 10;
        let tips: Vec<CladeFp> = (0..n).map(|_| CladeFp::random(&mut rng)).collect();

        let mut seen = HashSet::new();
        // Enumerate all 2^n - 1 non-empty subsets
        for mask in 1u32..(1 << n) {
            let mut fp = CladeFp::empty();
            for i in 0..n {
                if mask & (1 << i) != 0 {
                    fp = fp.union(&tips[i]);
                }
            }
            assert!(
                seen.insert(fp),
                "Collision: subset mask {:b} produced a duplicate fingerprint",
                mask
            );
        }
        assert_eq!(seen.len(), (1 << n) - 1);
    }

    #[test]
    fn assign_tip_fps_returns_distinct_nonempty_fps_for_each_tip() {
        use crate::newick::NewickNode;
        use crate::refs::{Pool, TestPool};

        let pool = TestPool::new();
        let tree = NewickTree::new(pool.alloc(NewickNode::inner_node(
            "root",
            0.0,
            vec![
                pool.alloc(NewickNode::inner_node(
                    "",
                    0.0,
                    vec![
                        pool.alloc(NewickNode::leaf("A", 0.0)),
                        pool.alloc(NewickNode::leaf("B", 0.0)),
                    ],
                )),
                pool.alloc(NewickNode::leaf("C", 0.0)),
            ],
        )));

        let mut rng = Pcg64Mcg::seed_from_u64(42);
        let tip_fps = assign_tip_fps(&tree, &mut rng);

        assert_eq!(tip_fps.len(), 3);
        assert!(tip_fps.contains_key("A"));
        assert!(tip_fps.contains_key("B"));
        assert!(tip_fps.contains_key("C"));

        let values: Vec<CladeFp> = tip_fps.values().copied().collect();
        for fp in &values {
            assert!(!fp.is_empty());
        }
        let unique: HashSet<CladeFp> = values.into_iter().collect();
        assert_eq!(unique.len(), 3);
    }

    #[test]
    fn rf_dist_identical_is_zero() {
        let a = vec![CladeFp::new(1), CladeFp::new(3), CladeFp::new(7)];
        assert_eq!(calc_rf_dist(&a, &a), 0);
    }

    #[test]
    fn rf_dist_completely_disjoint() {
        let a = vec![CladeFp::new(1), CladeFp::new(3)];
        let b = vec![CladeFp::new(2), CladeFp::new(4)];
        assert_eq!(calc_rf_dist(&a, &b), 4);
    }

    #[test]
    fn rf_dist_partial_overlap() {
        let a = vec![CladeFp::new(1), CladeFp::new(3), CladeFp::new(4), CladeFp::new(7)];
        let b = vec![CladeFp::new(2), CladeFp::new(3), CladeFp::new(5), CladeFp::new(7)];
        assert_eq!(calc_rf_dist(&a, &b), 4);
    }
}
