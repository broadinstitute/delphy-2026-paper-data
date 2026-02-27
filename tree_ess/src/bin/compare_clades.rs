//! Compare clades across two BEAST-like .trees files.
//!
//! For every clade (set of descendant tips) that appears in any tree in either file,
//! computes posterior support and mean MRCA date, with standard errors calibrated using
//! a tree topology ESS (Magee et al 2024).
//!
//! Inspired by BEAST 2's CladeSetComparator (https://www.beast2.org/2020/04/20/comparing-tree-sets.html).
//!
//! Output is JSON suitable for scatter-plotting support_A vs support_B and date_A vs date_B.

extern crate tree_ess;

use clap::Parser;
use itertools::Itertools;
use ndarray::{Array1, Array2};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io;
use tree_ess::burnin::BurninSpec;
use tree_ess::clade_fp::{assign_tip_fps, calc_rf_dist, CladeFp};
use tree_ess::newick::NewickTree;
use tree_ess::nexus_reader::NexusReader;
use tree_ess::refs::AllocPool;
use tree_ess::trees::{NodeLike, TraversalAction, TreeLike};

#[derive(Parser)]
#[command(
    version,
    about = "Compare clade support and dates across two .trees files"
)]
#[command(name = "compare_clades")]
struct Args {
    /// First .trees file
    file_a: String,

    /// Second .trees file
    file_b: String,

    /// Percentage of initial trees to discard as burn-in
    #[arg(long, group = "burn-in")]
    burnin_pct: Option<f64>,

    /// Number of initial trees to discard as burn-in
    #[arg(long, group = "burn-in")]
    burnin_trees: Option<usize>,

    /// Only output clades with support >= this threshold in at least one file (default: 0)
    #[arg(long, default_value = "0")]
    min_support: f64,
}

// -- Date parsing --

fn parse_tip_date(name: &str) -> Option<f64> {
    let date_str = name.rsplit('|').next()?;
    let parts: Vec<&str> = date_str.split('-').collect();
    if parts.len() != 3 {
        return None; // Only full YYYY-MM-DD dates
    }
    let year: i32 = parts[0].parse().ok()?;
    let month: u32 = parts[1].parse().ok()?;
    let day: u32 = parts[2].parse().ok()?;
    if !(1..=12).contains(&month) || !(1..=31).contains(&day) {
        return None;
    }
    Some(decimal_date(year, month, day))
}

fn decimal_date(year: i32, month: u32, day: u32) -> f64 {
    let is_leap = (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0);
    let days_in_year: f64 = if is_leap { 366.0 } else { 365.0 };
    let days_before_month: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let mut doy = days_before_month[(month - 1) as usize] + day;
    if is_leap && month > 2 {
        doy += 1;
    }
    year as f64 + (doy as f64 - 1.0) / days_in_year
}

// -- Clade definitions --

#[derive(Debug, Serialize, Clone, Eq, PartialEq)]
enum CladeDefinition {
    TipClade {
        name: String,
        fp: CladeFp,
    },
    InnerNodeClade {
        subclade1: CladeFp,
        subclade2: CladeFp,
        size: usize,
    },
}

impl CladeDefinition {
    fn fp(&self) -> CladeFp {
        match self {
            CladeDefinition::TipClade { fp, .. } => fp.clone(),
            CladeDefinition::InnerNodeClade {
                subclade1,
                subclade2,
                ..
            } => subclade1.union(&subclade2),
        }
    }

    fn size(&self) -> usize {
        match self {
            CladeDefinition::TipClade { .. } => 1,
            CladeDefinition::InnerNodeClade { size, .. } => *size,
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

type CladeMap = HashMap<CladeFp, CladeDefinition>;

fn tips_in_clade(fp: &CladeFp, clade_map: &CladeMap) -> Vec<String> {

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

// -- Per-tree processing --

struct CladeInTreeInfo {
    date: f64,
}

struct TreeResult {
    sorted_clade_fps: Vec<CladeFp>,
    clade_fps_2_info: HashMap<CladeFp, CladeInTreeInfo>,
}

fn find_root_date(tree: &NewickTree, exact_tip_dates: &HashMap<String, f64>) -> f64 {
    // Use the fixed date of some tip to back out the date of the root.
    // Panics if no tip has an exact date

    // We put up with the possibility of some round-off error in adding and subtracting
    // branch lengths as we traverse the tree.  The working assumption is that most tips
    // have exact dates, so it's very likely that just the initial descent from the root
    // to the first tip suffices to nail down the root date.  If that fails, then likely
    // one of the next few tips visited will have an exact date, thus producing a root date.
    // Hence, the amount of roundoff that we expect in practice is very small.

    let mut root_to_node_dist = 0.0;

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                if node_ref != tree.root {
                    // root-to-root distance is 0 by definition; ignore branch length then
                    root_to_node_dist += node.branch_length
                };

                if node.is_tip()
                    && let Some(exact_tip_date) = exact_tip_dates.get(node.name.as_str())
                {
                    // Found a tip with an exact date: done!
                    let root_date = exact_tip_date - root_to_node_dist;
                    return root_date;
                }
            }
            TraversalAction::Exit => {
                // If we get here when `node` is the root, then we're going to fail
                // anyway, so there's nothing to gain from subtracting `node.branch_length`
                // only for non-root nodes.
                root_to_node_dist -= node.branch_length;
            }
        }
    }

    panic!("No tip found with an exact date!  Cannot deduce root date!");
}

fn process_tree(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    exact_tip_dates: &HashMap<String, f64>,
    clade_map: &mut CladeMap,
) -> TreeResult {

    let root_date = find_root_date(tree, exact_tip_dates);

    let mut node_clade_defn = CladeDefinition::InnerNodeClade {
        subclade1: CladeFp::empty(),
        subclade2: CladeFp::empty(),
        size: 0,
    };
    let mut node_clade_defn_stack: Vec<CladeDefinition> = vec![];

    let mut node_date = root_date;
    let mut node_date_stack: Vec<f64> = vec![];

    let mut clade_fps_2_info: HashMap<CladeFp, CladeInTreeInfo> = HashMap::new();

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                // Shift focus from parent to node

                node_date_stack.push(node_date);
                if node_ref != tree.root {
                    // root date is root_date by definition; ignore branch length then
                    node_date += node.branch_length
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
                            size: 0,
                        }
                    }
                };
                node_clade_defn_stack.push(std::mem::replace(&mut node_clade_defn, new_clade_defn));

                // Invariant: at this point, node_date and node_clade_defn refer to the
                // date and (partial) fingerprint for the current node after entering it
            }
            TraversalAction::Exit => {
                // Invariant: at this point, node_date and node_clade_defn refer to the
                // date and (partial) fingerprint for the current node before exiting it

                assert!(node_clade_defn.is_complete());
                let node_fp = node_clade_defn.fp();

                clade_fps_2_info.insert(node_fp, CladeInTreeInfo { date: node_date });

                // After this line, node_clade_defn becomes the parent's clade definition
                let child_clade_defn =
                    std::mem::replace(&mut node_clade_defn, node_clade_defn_stack.pop().unwrap());
                let child_clade_size = child_clade_defn.size();
                clade_map.entry(node_fp).or_insert(child_clade_defn); // Save one example decomposition of this clade for later reconstruction

                // Shift focus from node to parent, merging this node's fingerprint into it

                // node_clade_defn is already the parent's clade definition; merge exited child's
                // fingerprint if necessary
                node_clade_defn = match node_clade_defn {
                    CladeDefinition::TipClade { .. } => node_clade_defn,
                    CladeDefinition::InnerNodeClade {
                        subclade1,
                        subclade2,
                        size,
                    } => {
                        if subclade1.is_empty() {
                            CladeDefinition::InnerNodeClade {
                                subclade1: node_fp,
                                subclade2: CladeFp::empty(),
                                size: child_clade_size,
                            }
                        } else if subclade2.is_empty() {
                            CladeDefinition::InnerNodeClade {
                                subclade1,
                                subclade2: node_fp,
                                size: size + child_clade_size,
                            }
                        } else {
                            panic!("Non-bifurcating tree?")
                        }
                    }
                };

                let parent_date = node_date_stack.pop().unwrap();
                node_date = parent_date;
            }
        }
    }

    let sorted_clade_fps: Vec<CladeFp> = clade_fps_2_info.keys().cloned().sorted().collect();

    TreeResult {
        sorted_clade_fps,
        clade_fps_2_info,
    }
}

// -- Clade accumulator --

#[derive(Default)]
struct CladeAccumulator {
    count: usize,
    date_sum: f64,
    date_sq_sum: f64,
}

// -- ESS computation (adapted from calc_tree_ess.rs) --

fn compute_ess(sorted_clade_fps: &[Vec<CladeFp>]) -> f64 {
    let n = sorted_clade_fps.len();
    if n <= 1 {
        return n as f64;
    }

    // All-pairs RF distance matrix
    let mut rf_dist = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        eprintln!("  Computing RF distances from tree {} of {}...", i + 1, n);
        for j in (i + 1)..n {
            let d =
                calc_rf_dist(&sorted_clade_fps[i], &sorted_clade_fps[j]) as f64;
            rf_dist[(i, j)] = d;
            rf_dist[(j, i)] = d;
        }
    }

    // Var[tau] = 1/(n(n-1)) sum_{i<j} d^2
    let mut var_tau = 0.0;
    for &d in rf_dist.iter() {
        var_tau += d * d;
    }
    var_tau /= 2.0 * (n as f64) * ((n - 1) as f64);

    if var_tau == 0.0 {
        return n as f64; // All trees identical
    }

    // E[Delta_s^2] at each lag s
    let mut mean_delta_s_2 = Array1::<f64>::zeros(n);
    for s in 0..n {
        for i in 0..(n - s) {
            let d = rf_dist[(i, i + s)];
            mean_delta_s_2[s] += d * d;
        }
        mean_delta_s_2[s] /= (n - s) as f64;
    }

    let rho_s = 1.0 - 0.5 * mean_delta_s_2 / var_tau;

    // ACT = 1 + 2 * sum_{s=1}^{smax-1} rho_s[s]
    // smax = first s where rho_s[s-1] + rho_s[s] <= 0
    let mut smax = 1;
    while (smax < n) && ((rho_s[smax - 1] + rho_s[smax]) > 0.0) {
        smax += 1;
    }
    let mut act = 1.0;
    for s in 1..smax {
        act += 2.0 * rho_s[s];
    }

    n as f64 / act
}

// -- Output types --

#[derive(Serialize)]
struct FileInfo {
    path: String,
    num_trees_total: usize,
    num_trees_burnin: usize,
    num_trees_post_burnin: usize,
    tree_ess: f64,
}

#[derive(Serialize)]
struct CladeFileStats {
    support: f64,
    support_se: f64,
    mean_date: f64,
    date_se: f64,
}

#[derive(Serialize)]
struct CladeResult {
    clade_fp: CladeFp,
    clade_definition: CladeDefinition,
    tips_in_clade: Vec<String>,
    file_a: CladeFileStats,
    file_b: CladeFileStats,
}

#[derive(Serialize)]
struct Output {
    file_a: FileInfo,
    file_b: FileInfo,
    num_tips: usize,
    clades: Vec<CladeResult>,
}

// -- File processing --

struct FileProcessingResult {
    accumulators: HashMap<CladeFp, CladeAccumulator>,
    ess: f64,
    num_total: usize,
    num_burnin: usize,
    num_post_burnin: usize,
}

fn process_trees(
    file_path: &str,
    trees: &Vec<(u64, NewickTree)>,
    burnin: &BurninSpec,
    tip_fps: &HashMap<String, CladeFp>,
    exact_tip_dates: &HashMap<String, f64>,
    clade_map: &mut CladeMap,
) -> Result<FileProcessingResult, Box<dyn Error>> {
    let num_total = trees.len();
    assert!(num_total >= 1, "No trees in {}", file_path);

    let num_burnin = burnin.first_sample_idx(num_total);
    let num_post_burnin = num_total - num_burnin;

    eprintln!(
        "  {} trees total, {} burn-in, {} post-burn-in",
        num_total, num_burnin, num_post_burnin
    );

    let mut accumulators: HashMap<CladeFp, CladeAccumulator> = HashMap::new();
    let mut all_sorted_clade_fps: Vec<Vec<CladeFp>> = Vec::new();

    for (i, (_state, tree)) in trees.into_iter().enumerate() {
        if i < num_burnin {
            continue;
        }

        let result = process_tree(&tree, tip_fps, exact_tip_dates, clade_map);

        // Accumulate per-clade statistics
        for (&clade_fp, info) in &result.clade_fps_2_info {
            let acc = accumulators.entry(clade_fp).or_default();
            acc.count += 1;
            acc.date_sum += info.date;
            acc.date_sq_sum += info.date * info.date;
        }

        all_sorted_clade_fps.push(result.sorted_clade_fps);
    }

    eprintln!("Computing tree ESS for {}...", file_path);
    let ess = compute_ess(&all_sorted_clade_fps);
    eprintln!("  ESS = {:.1}", ess);

    Ok(FileProcessingResult {
        accumulators,
        ess,
        num_total,
        num_burnin,
        num_post_burnin,
    })
}

fn make_clade_file_stats(
    acc: Option<&CladeAccumulator>,
    num_post_burnin: usize,
    ess: f64,
) -> CladeFileStats {
    match acc {
        None | Some(CladeAccumulator { count: 0, .. }) => CladeFileStats {
            support: 0.0,
            support_se: 0.0,
            mean_date: f64::NAN,
            date_se: f64::NAN,
        },
        Some(acc) => {
            let p = acc.count as f64 / num_post_burnin as f64;
            let support_se = (p * (1.0 - p) / ess).sqrt();

            let mean_date = acc.date_sum / acc.count as f64;
            let variance = (acc.date_sq_sum / acc.count as f64) - mean_date * mean_date;
            let sd = variance.max(0.0).sqrt();
            let ess_for_height = ess.min(acc.count as f64);
            let date_se = sd / ess_for_height.sqrt();

            CladeFileStats {
                support: p,
                support_se,
                mean_date,
                date_se,
            }
        }
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    let pool = AllocPool::new();
    let args = Args::parse();

    let burnin_spec = BurninSpec::from_options(args.burnin_pct, args.burnin_trees);

    // Since the clade fingerprints make it into output files, we use a PRNG with a fixed seed
    // (doesn't affect the chances of collisions).  We use Pcg64Mcg specifically because it is
    // guaranteed to be portable, i.e., produces the same stream of numbers on every platform and
    // on every future version of the `rand` family of crates.
    // See https://docs.rs/rand_pcg/latest/rand_pcg/ for details
    let mut rng = Pcg64Mcg::seed_from_u64(1);

    // Parse the input trees from both files
    eprintln!("Reading trees in {}...", args.file_a);
    let file = File::open(&args.file_a)?;
    let reader = BufReader::new(file);
    let file_a_trees = NexusReader::parse(reader, 0, 1, &pool)?;
    assert!(!file_a_trees.is_empty(), "No trees in {}", args.file_a);

    eprintln!("Reading trees in {}...", args.file_b);
    let file = File::open(&args.file_b)?;
    let reader = BufReader::new(file);
    let file_b_trees = NexusReader::parse(reader, 0, 1, &pool)?;
    assert!(!file_b_trees.is_empty(), "No trees in {}", args.file_b);

    // Extract tip names from first tree of file A, then check that all the other trees
    // have the same tip names
    eprintln!("Extracting tip names from first tree in {}...", args.file_a);
    let (first_state_a, first_tree_a) = file_a_trees.first().unwrap();
    let tip_fps = assign_tip_fps(first_tree_a, &mut rng);
    let num_tips = tip_fps.len();
    eprintln!("  {} tips found", num_tips);

    eprintln!("Checking that all trees have matching tips...");
    for (trees, name) in [(&file_a_trees, &args.file_a), (&file_b_trees, &args.file_b)].iter() {
        for (state, tree) in trees.iter() {

            let this_trees_tips: HashSet<String> =
                tree.any_order_iter()
                    .filter(|node_ref| node_ref.borrow().is_tip())
                    .map(|node_ref| node_ref.borrow().name.clone())
                    .collect();

            let tips_are_equal =
                this_trees_tips.len() == num_tips &&
                    this_trees_tips.iter().all(|tip_name| tip_fps.contains_key(tip_name));

            if !tips_are_equal {
                panic!(
                    "Tip mismatch between state {} in tree in {} and state {} in {}",
                    first_state_a, args.file_a, state, name
                );
            }
        }
    }

    // Build tip date map
    let mut exact_tip_dates: HashMap<String, f64> = HashMap::new();
    for name in tip_fps.keys() {
        if let Some(date) = parse_tip_date(name) {
            exact_tip_dates.insert(name.clone(), date);
        }
    }
    eprintln!("  {} tips with full dates", exact_tip_dates.len());
    if exact_tip_dates.is_empty() {
        panic!("WARNING: No tips with full YYYY-MM-DD dates found.");
    }

    // Clade definitions will accumulate here
    let mut clade_map: CladeMap = HashMap::new();

    // Process both files
    let result_a = process_trees(
        &args.file_a,
        &file_a_trees,
        &burnin_spec,
        &tip_fps,
        &exact_tip_dates,
        &mut clade_map,
    )?;
    let result_b = process_trees(
        &args.file_b,
        &file_b_trees,
        &burnin_spec,
        &tip_fps,
        &exact_tip_dates,
        &mut clade_map,
    )?;

    // Collect all clade fingerprints
    let mut all_fps: Vec<CladeFp> = clade_map.keys().copied().collect();
    all_fps.sort();

    // Build output, applying min-support filter
    eprintln!(
        "Building output ({} clades total, min-support={})...",
        all_fps.len(),
        args.min_support
    );
    let mut clades: Vec<CladeResult> = Vec::new();
    for fp in &all_fps {
        let stats_a = make_clade_file_stats(
            result_a.accumulators.get(fp),
            result_a.num_post_burnin,
            result_a.ess,
        );
        let stats_b = make_clade_file_stats(
            result_b.accumulators.get(fp),
            result_b.num_post_burnin,
            result_b.ess,
        );

        if stats_a.support < args.min_support && stats_b.support < args.min_support {
            continue;
        }

        clades.push(CladeResult {
            clade_fp: fp.clone(),
            clade_definition: clade_map.get(fp).unwrap().clone(),
            tips_in_clade: tips_in_clade(&fp, &clade_map),
            file_a: stats_a,
            file_b: stats_b,
        });
    }
    eprintln!("  {} clades pass min-support filter", clades.len());

    let output = Output {
        file_a: FileInfo {
            path: args.file_a,
            num_trees_total: result_a.num_total,
            num_trees_burnin: result_a.num_burnin,
            num_trees_post_burnin: result_a.num_post_burnin,
            tree_ess: result_a.ess,
        },
        file_b: FileInfo {
            path: args.file_b,
            num_trees_total: result_b.num_total,
            num_trees_burnin: result_b.num_burnin,
            num_trees_post_burnin: result_b.num_post_burnin,
            tree_ess: result_b.ess,
        },
        num_tips,
        clades,
    };

    serde_json::to_writer_pretty(io::stdout(), &output)?;
    eprintln!("Done!");

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::CladeDefinition::{InnerNodeClade, TipClade};
    use itertools::Itertools;
    use tree_ess::newick::NewickNode;
    use tree_ess::refs::Pool;

    #[test]
    fn test_parse_tip_date_full() {
        assert!((parse_tip_date("EBOV-20547|2015-09-07").unwrap() - 2015.6822).abs() < 0.01);
    }

    #[test]
    fn test_parse_tip_date_partial_month() {
        assert_eq!(parse_tip_date("OP415257|9000312|UK|2022-08"), None);
    }

    #[test]
    fn test_parse_tip_date_no_pipe() {
        assert_eq!(parse_tip_date("just_a_name"), None);
    }

    #[test]
    fn test_decimal_date() {
        // Jan 1 of a non-leap year
        assert!((decimal_date(2023, 1, 1) - 2023.0).abs() < 1e-6);
        // Dec 31 of a non-leap year
        assert!((decimal_date(2023, 12, 31) - (2023.0 + 364.0 / 365.0)).abs() < 1e-6);
    }

    #[test]
    fn test_decimal_date_leap_year() {
        // March 1 of a leap year
        assert!((decimal_date(2024, 3, 1) - (2024.0 + 60.0 / 366.0)).abs() < 1e-6);
    }

    #[test]
    fn test_process_tree_basic() {
        use tree_ess::refs::TestPool;

        let (name_a, fp_a) = ("A|2020-01-01", CladeFp::new(0b001));
        let (name_b, fp_b) = ("B|2020-07-01", CladeFp::new(0b010));
        let (name_c, fp_c) = ("C|2021-01-01", CladeFp::new(0b100));

        let pool = TestPool::new();
        let tree = NewickTree::new(pool.alloc(NewickNode::inner_node(
            "root",
            0.0,
            vec![
                pool.alloc(NewickNode::inner_node(
                    "",
                    1.0,
                    vec![
                        pool.alloc(NewickNode::leaf(name_a, 0.5)),
                        pool.alloc(NewickNode::leaf(name_b, 0.5)),
                    ],
                )),
                pool.alloc(NewickNode::leaf(name_c, 2.0)),
            ],
        )));

        let tip_fps: HashMap<String, CladeFp> = HashMap::from([
            (name_a.to_string(), fp_a),
            (name_b.to_string(), fp_b),
            (name_c.to_string(), fp_c),
        ]);
        let exact_tip_dates: HashMap<String, f64> = HashMap::from([
            (name_a.to_string(), parse_tip_date(name_a).unwrap()),
            (name_b.to_string(), parse_tip_date(name_b).unwrap()),
            (name_c.to_string(), parse_tip_date(name_c).unwrap()),
        ]);

        let mut clade_map: CladeMap = HashMap::new();

        let result = process_tree(&tree, &tip_fps, &exact_tip_dates, &mut clade_map);

        // 2 inner nodes: {A,B} and the root {A,B,C}
        let fp_ab = fp_a.union(&fp_b);
        let fp_abc = fp_ab.union(&fp_c);
        assert_eq!(
            result.sorted_clade_fps,
            vec![fp_a, fp_b, fp_c, fp_ab, fp_abc]
                .iter()
                .copied()
                .sorted()
                .collect::<Vec<CladeFp>>()
        );

        // Should have discovered clades
        assert_eq!(
            clade_map,
            HashMap::from([
                (
                    fp_a,
                    TipClade {
                        name: name_a.to_string(),
                        fp: fp_a
                    }
                ),
                (
                    fp_b,
                    TipClade {
                        name: name_b.to_string(),
                        fp: fp_b
                    }
                ),
                (
                    fp_c,
                    TipClade {
                        name: name_c.to_string(),
                        fp: fp_c
                    }
                ),
                (
                    fp_ab,
                    InnerNodeClade {
                        subclade1: fp_a,
                        subclade2: fp_b,
                        size: 2
                    }
                ),
                (
                    fp_abc,
                    InnerNodeClade {
                        subclade1: fp_ab,
                        subclade2: fp_c,
                        size: 3
                    }
                ),
            ])
        );

        // Reconstruct tip names from clade fingerprints
        assert_eq!(tips_in_clade(&fp_a, &clade_map), vec![name_a.to_string()]);
        assert_eq!(tips_in_clade(&fp_b, &clade_map), vec![name_b.to_string()]);
        assert_eq!(tips_in_clade(&fp_c, &clade_map), vec![name_c.to_string()]);
        assert_eq!(
            tips_in_clade(&fp_ab, &clade_map),
            vec![name_a.to_string(), name_b.to_string()]
        );
        assert_eq!(
            tips_in_clade(&fp_abc, &clade_map),
            vec![name_a.to_string(), name_b.to_string(), name_c.to_string()]
        );
    }
}
