//! Calculate clade coverage statistics for one replicate of a well-calibrated simulation study.
//!
//! Given a true (generating) tree and a set of posterior trees from inference, bins posterior
//! clades by their posterior support and reports how many are actually present in the true tree.
//!
//! Designed to be called once per replicate, outputting one TSV line per call.  A downstream
//! script aggregates across replicates to build the clade coverage bar chart.
//!
//! Inspired by BEAST 2's CladeCoverageCalculator:
//! https://github.com/christiaanjs/beast-validation/blob/master/src/beastvalidation/experimenter/CladeCoverageCalculator.java

extern crate tree_ess;

use clap::Parser;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::cmp;
use tree_ess::newick::NewickTree;
use tree_ess::nexus_reader::NexusReader;
use tree_ess::refs::AllocPool;
use tree_ess::trees::{NodeLike, TraversalAction, TreeLike};

#[derive(Parser)]
#[command(
    version,
    about = "Calculate clade coverage statistics for one WCSS replicate"
)]
#[command(name = "calc_clade_coverage")]
struct Args {
    /// NEXUS file containing the true (generating) tree
    true_tree: String,

    /// BEAST/Delphy .trees file with posterior tree samples
    posterior_trees: String,

    /// Replicate identifier (included as first column of output)
    #[arg(long)]
    replicate: String,

    /// Percentage of posterior trees to discard as burn-in
    #[arg(long, group = "burn-in")]
    burnin_pct: Option<f64>,

    /// Number of posterior trees to discard as burn-in
    #[arg(long, group = "burn-in")]
    burnin_trees: Option<usize>,

    /// Number of equal-width bins spanning [0, 1] (default: 10 = 10% bins)
    #[arg(long, default_value = "10")]
    num_bins: usize,

    /// Prepend a TSV header line before the data line
    #[arg(long)]
    header: bool,
}

// -- Burn-in (same pattern as calc_tree_ess.rs and compare_clades.rs) --

#[derive(Debug)]
enum BurninSpec {
    Fract(f64),
    Trees(usize),
}
impl BurninSpec {
    fn first_sample_idx(&self, num_trees: usize) -> usize {
        match *self {
            BurninSpec::Fract(pct) => {
                assert!((0.0..=1.0).contains(&pct));
                (num_trees as f64 * pct).floor() as usize
            }
            BurninSpec::Trees(burnin_trees) => cmp::min(num_trees, burnin_trees),
        }
    }
}

// -- Clade fingerprint (same as compare_clades.rs) --

#[derive(Debug, Eq, PartialEq, Ord, PartialOrd, Hash, Clone, Copy)]
struct CladeFp(u64);
impl CladeFp {
    fn empty() -> CladeFp {
        CladeFp(0)
    }

    fn random(rng: &mut dyn Rng) -> CladeFp {
        CladeFp(rng.next_u64())
    }

    fn union(&self, other: &CladeFp) -> CladeFp {
        CladeFp(self.0 ^ other.0)
    }
}

// -- Tip fingerprint assignment (same pattern as compare_clades.rs) --

fn assign_tip_fps(tree: &NewickTree, rng: &mut dyn Rng) -> HashMap<String, CladeFp> {
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

/// Traverse a tree and return clade fingerprints for all non-trivial clades
/// (excluding tip singletons and the root clade).
fn extract_nontrivial_clades(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
) -> Vec<CladeFp> {
    let mut fp_stack: Vec<CladeFp> = Vec::new();
    let mut result: Vec<CladeFp> = Vec::new();

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                if node.is_tip() {
                    let fp = *tip_fps
                        .get(node.name.as_str())
                        .expect("Unknown tip name");
                    fp_stack.push(fp);
                } else {
                    fp_stack.push(CladeFp::empty());
                }
            }
            TraversalAction::Exit => {
                let fp = fp_stack.pop().expect("Improper nesting?");

                // Exclude tips and the root
                let is_tip = node.is_tip();
                let is_root = node_ref == tree.root;
                if !is_tip && !is_root {
                    result.push(fp);
                }

                if let Some(parent_fp) = fp_stack.last_mut() {
                    *parent_fp = parent_fp.union(&fp);
                }
            }
        }
    }

    result
}

// -- Binning --

fn compute_bin_index(support: f64, bin_width: f64, num_bins: usize) -> usize {
    let b = (support / bin_width).floor() as usize;
    cmp::min(b, num_bins - 1)
}

fn main() -> Result<(), Box<dyn Error>> {
    let pool = AllocPool::new();
    let args = Args::parse();

    let burnin_spec = if let Some(pct) = args.burnin_pct {
        assert!((0.0..=100.0).contains(&pct));
        BurninSpec::Fract(pct / 100.0)
    } else if let Some(burnin_trees) = args.burnin_trees {
        BurninSpec::Trees(burnin_trees)
    } else {
        BurninSpec::Fract(0.10)
    };

    let num_bins = args.num_bins;
    assert!(num_bins >= 1, "num-bins must be at least 1");
    let bin_width = 1.0 / num_bins as f64;

    // Deterministic fingerprints (same pattern as compare-clades)
    let mut rng = Pcg64Mcg::seed_from_u64(1);

    // Step 1: Parse inputs

    eprintln!("Reading true tree from {}...", args.true_tree);
    let file = File::open(&args.true_tree)?;
    let reader = BufReader::new(file);
    let true_trees = NexusReader::parse(reader, 0, 1, &pool)?;
    assert!(
        true_trees.len() == 1,
        "Expected exactly 1 tree in {}, found {}",
        args.true_tree,
        true_trees.len()
    );
    let true_tree = &true_trees[0].1;

    eprintln!("Reading posterior trees from {}...", args.posterior_trees);
    let file = File::open(&args.posterior_trees)?;
    let reader = BufReader::new(file);
    let posterior_trees = NexusReader::parse(reader, 0, 1, &pool)?;
    assert!(
        !posterior_trees.is_empty(),
        "No trees in {}",
        args.posterior_trees
    );

    // Assign fingerprints from the first posterior tree's tips
    let tip_fps = assign_tip_fps(&posterior_trees[0].1, &mut rng);
    let num_tips = tip_fps.len();
    eprintln!("  {} tips found", num_tips);

    // Verify the true tree has the same tips
    let true_tree_tips: HashSet<String> = true_tree
        .any_order_iter()
        .filter(|node_ref| node_ref.borrow().is_tip())
        .map(|node_ref| node_ref.borrow().name.clone())
        .collect();
    assert!(
        true_tree_tips.len() == num_tips
            && true_tree_tips.iter().all(|name| tip_fps.contains_key(name)),
        "Tip mismatch between true tree and posterior trees"
    );

    // Apply burn-in
    let num_total = posterior_trees.len();
    let num_burnin = burnin_spec.first_sample_idx(num_total);
    let num_post_burnin = num_total - num_burnin;
    eprintln!(
        "  {} trees total, {} burn-in, {} post-burn-in",
        num_total, num_burnin, num_post_burnin
    );
    assert!(num_post_burnin > 0, "No trees remaining after burn-in");

    // Step 2: Compute posterior clade support
    eprintln!("Computing posterior clade support...");
    let mut clade_counts: HashMap<CladeFp, usize> = HashMap::new();

    for (i, (_state, tree)) in posterior_trees.iter().enumerate() {
        if i < num_burnin {
            continue;
        }
        let clades = extract_nontrivial_clades(tree, &tip_fps);
        for &fp in &clades {
            *clade_counts.entry(fp).or_insert(0) += 1;
        }
    }
    eprintln!(
        "  {} distinct non-trivial clades in posterior",
        clade_counts.len()
    );

    // Step 3: Merge true tree clades into posterior support map
    eprintln!("Extracting true tree clades...");
    let true_clade_fps: HashSet<CladeFp> =
        extract_nontrivial_clades(true_tree, &tip_fps).into_iter().collect();
    eprintln!("  {} non-trivial clades in true tree", true_clade_fps.len());

    // Ensure every true clade is in the counts map (with 0 if not seen in posterior)
    for &fp in &true_clade_fps {
        clade_counts.entry(fp).or_insert(0);
    }

    // Step 4: Bin and count
    let mut totals = vec![0usize; num_bins];
    let mut true_hits = vec![0usize; num_bins];

    // Totals: iterate over all clades in the HashMap
    for &count in clade_counts.values() {
        let support = count as f64 / num_post_burnin as f64;
        let b = compute_bin_index(support, bin_width, num_bins);
        totals[b] += 1;
    }

    // True hits: iterate over true clades only
    for &fp in &true_clade_fps {
        let count = *clade_counts.get(&fp).unwrap();
        let support = count as f64 / num_post_burnin as f64;
        let b = compute_bin_index(support, bin_width, num_bins);
        true_hits[b] += 1;
    }

    // Step 5: Output
    let bin_boundaries: Vec<(usize, usize)> = (0..num_bins)
        .map(|b| {
            let lo = (b as f64 * bin_width * 100.0).round() as usize;
            let hi = ((b + 1) as f64 * bin_width * 100.0).round() as usize;
            (lo, hi)
        })
        .collect();

    if args.header {
        let mut header_parts = vec!["replicate".to_string()];
        for &(lo, hi) in &bin_boundaries {
            header_parts.push(format!("totals_{}_{}", lo, hi));
            header_parts.push(format!("true_hits_{}_{}", lo, hi));
            header_parts.push(format!("frac_{}_{}", lo, hi));
        }
        println!("{}", header_parts.join("\t"));
    }

    let mut data_parts = vec![args.replicate.clone()];
    for b in 0..num_bins {
        data_parts.push(format!("{}", totals[b]));
        data_parts.push(format!("{}", true_hits[b]));
        if totals[b] > 0 {
            data_parts.push(format!("{:.6}", true_hits[b] as f64 / totals[b] as f64));
        } else {
            data_parts.push("NaN".to_string());
        }
    }
    println!("{}", data_parts.join("\t"));

    eprintln!("Done!");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use tree_ess::newick::NewickNode;
    use tree_ess::refs::{Pool, TestPool};

    fn make_tip_fps() -> HashMap<String, CladeFp> {
        HashMap::from([
            ("A".to_string(), CladeFp(0b0001)),
            ("B".to_string(), CladeFp(0b0010)),
            ("C".to_string(), CladeFp(0b0100)),
            ("D".to_string(), CladeFp(0b1000)),
        ])
    }

    // ((A,B),(C,D)) — clades: {A,B} and {C,D}
    fn make_tree_ab_cd(pool: &dyn Pool<NewickNode>) -> NewickTree {
        NewickTree::new(pool.alloc(NewickNode::inner_node(
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
                pool.alloc(NewickNode::inner_node(
                    "",
                    0.0,
                    vec![
                        pool.alloc(NewickNode::leaf("C", 0.0)),
                        pool.alloc(NewickNode::leaf("D", 0.0)),
                    ],
                )),
            ],
        )))
    }

    // ((A,C),(B,D)) — clades: {A,C} and {B,D}
    fn make_tree_ac_bd(pool: &dyn Pool<NewickNode>) -> NewickTree {
        NewickTree::new(pool.alloc(NewickNode::inner_node(
            "root",
            0.0,
            vec![
                pool.alloc(NewickNode::inner_node(
                    "",
                    0.0,
                    vec![
                        pool.alloc(NewickNode::leaf("A", 0.0)),
                        pool.alloc(NewickNode::leaf("C", 0.0)),
                    ],
                )),
                pool.alloc(NewickNode::inner_node(
                    "",
                    0.0,
                    vec![
                        pool.alloc(NewickNode::leaf("B", 0.0)),
                        pool.alloc(NewickNode::leaf("D", 0.0)),
                    ],
                )),
            ],
        )))
    }

    #[test]
    fn test_identical_trees() {
        let pool = TestPool::new();
        let tip_fps = make_tip_fps();

        let tree = make_tree_ab_cd(&pool);

        // Simulate 1 posterior tree identical to the true tree
        let mut clade_counts: HashMap<CladeFp, usize> = HashMap::new();
        for &fp in &extract_nontrivial_clades(&tree, &tip_fps) {
            *clade_counts.entry(fp).or_insert(0) += 1;
        }

        // True tree clades
        let true_clade_fps: HashSet<CladeFp> =
            extract_nontrivial_clades(&tree, &tip_fps).into_iter().collect();

        // Merge
        for &fp in &true_clade_fps {
            clade_counts.entry(fp).or_insert(0);
        }

        // Bin (10 bins, width 0.1)
        let num_bins = 10;
        let bin_width = 0.1;
        let num_post_burnin = 1;
        let mut totals = vec![0usize; num_bins];
        let mut true_hits = vec![0usize; num_bins];

        for &count in clade_counts.values() {
            let support = count as f64 / num_post_burnin as f64;
            let b = compute_bin_index(support, bin_width, num_bins);
            totals[b] += 1;
        }
        for &fp in &true_clade_fps {
            let count = *clade_counts.get(&fp).unwrap();
            let support = count as f64 / num_post_burnin as f64;
            let b = compute_bin_index(support, bin_width, num_bins);
            true_hits[b] += 1;
        }

        // All clades have support 1.0, so they land in the last bin
        assert_eq!(totals[9], 2); // {A,B} and {C,D}
        assert_eq!(true_hits[9], 2);
        for b in 0..9 {
            assert_eq!(totals[b], 0);
            assert_eq!(true_hits[b], 0);
        }
    }

    #[test]
    fn test_completely_different_trees() {
        let pool = TestPool::new();
        let tip_fps = make_tip_fps();

        let posterior_tree = make_tree_ab_cd(&pool);
        let true_tree = make_tree_ac_bd(&pool);

        // Posterior clades: {A,B} and {C,D}
        let mut clade_counts: HashMap<CladeFp, usize> = HashMap::new();
        for &fp in &extract_nontrivial_clades(&posterior_tree, &tip_fps) {
            *clade_counts.entry(fp).or_insert(0) += 1;
        }

        // True tree clades: {A,C} and {B,D}
        let true_clade_fps: HashSet<CladeFp> =
            extract_nontrivial_clades(&true_tree, &tip_fps).into_iter().collect();

        // No overlap
        for &fp in &true_clade_fps {
            assert!(!clade_counts.contains_key(&fp));
        }

        // Merge: true clades get count 0
        for &fp in &true_clade_fps {
            clade_counts.entry(fp).or_insert(0);
        }

        let num_bins = 10;
        let bin_width = 0.1;
        let num_post_burnin = 1;
        let mut totals = vec![0usize; num_bins];
        let mut true_hits = vec![0usize; num_bins];

        for &count in clade_counts.values() {
            let support = count as f64 / num_post_burnin as f64;
            let b = compute_bin_index(support, bin_width, num_bins);
            totals[b] += 1;
        }
        for &fp in &true_clade_fps {
            let count = *clade_counts.get(&fp).unwrap();
            let support = count as f64 / num_post_burnin as f64;
            let b = compute_bin_index(support, bin_width, num_bins);
            true_hits[b] += 1;
        }

        // Posterior clades {A,B} and {C,D} have support 1.0 → bin 9
        // True clades {A,C} and {B,D} have support 0.0 → bin 0
        assert_eq!(totals[0], 2);   // true-tree-only clades with 0 support
        assert_eq!(totals[9], 2);   // posterior clades with support 1.0
        assert_eq!(true_hits[0], 2); // true clades landed in bin 0
        assert_eq!(true_hits[9], 0); // posterior clades in bin 9 are not true
    }

    #[test]
    fn test_bin_boundaries() {
        let bin_width = 0.1;
        let num_bins = 10;

        // Support exactly 0.0 → bin 0
        assert_eq!(compute_bin_index(0.0, bin_width, num_bins), 0);

        // Support exactly 0.1 → bin 1 (not bin 0)
        assert_eq!(compute_bin_index(0.1, bin_width, num_bins), 1);

        // Support exactly 0.9 → bin 9
        assert_eq!(compute_bin_index(0.9, bin_width, num_bins), 9);

        // Support exactly 1.0 → capped to bin 9
        assert_eq!(compute_bin_index(1.0, bin_width, num_bins), 9);

        // Support 0.099999 → bin 0
        assert_eq!(compute_bin_index(0.099999, bin_width, num_bins), 0);

        // Support 0.5 → bin 5
        assert_eq!(compute_bin_index(0.5, bin_width, num_bins), 5);
    }
}
