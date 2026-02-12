//! Compare clades across two BEAST-like .trees files.
//!
//! For every clade (set of descendant tips) that appears in any tree in either file,
//! computes posterior support and mean MRCA date, with standard errors derived from
//! the global tree ESS (Magee et al 2024).
//!
//! Inspired by BEAST 2's CladeSetComparator (https://www.beast2.org/2020/04/20/comparing-tree-sets.html).
//!
//! Output is JSON suitable for scatter-plotting support_A vs support_B and date_A vs date_B.

extern crate tree_ess;

use clap::Parser;
use itertools::EitherOrBoth;
use ndarray::{Array1, Array2};
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::{cmp, io};
use tree_ess::newick::{NewickNode, NewickTree};
use tree_ess::nexus_reader::NexusReader;
use tree_ess::refs::{AllocPool, Pool};
use tree_ess::trees::{NodeLike, TraversalAction, TreeLike};

#[derive(Parser)]
#[command(version, about = "Compare clade support and dates across two .trees files")]
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

// -- Burn-in (adapted from calc_tree_ess.rs) --

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

// -- Tip fingerprints (from calc_tree_ess.rs) --

fn assign_tip_fingerprints(tree: &NewickTree) -> HashMap<String, u64> {
    let mut result = HashMap::new();
    for node_ref in tree.any_order_iter() {
        let node = node_ref.borrow();
        if node.is_tip() {
            assert!(!node.name.is_empty(), "Tip with no name!");
            result.insert(node.name.clone(), rand::random());
        }
    }
    result
}

// -- RF distance (from calc_tree_ess.rs) --

fn calc_rf_dist(
    sorted_split_fingerprints_a: &Vec<u64>,
    sorted_split_fingerprints_b: &Vec<u64>,
) -> u64 {
    use EitherOrBoth::{Both, Left, Right};
    let mut distance = 0;
    for pair in itertools::merge_join_by(
        sorted_split_fingerprints_a.iter(),
        sorted_split_fingerprints_b.iter(),
        |a, b| a.cmp(b),
    ) {
        match pair {
            Both(_, _) => (),
            Left(_) | Right(_) => distance += 1,
        }
    }
    distance
}

// -- Per-tree processing --

struct InnerNodeInfo {
    fingerprint: u64,
    dist_from_root: f64,
}

struct TreeResult {
    sorted_split_fingerprints: Vec<u64>,
    inner_nodes: Vec<InnerNodeInfo>,
    root_date: Option<f64>,
}

fn process_tree(
    tree: &NewickTree,
    tip_fingerprints: &HashMap<String, u64>,
    tip_dates: &HashMap<String, Option<f64>>,
    tip_name_to_index: &HashMap<String, usize>,
    clade_tip_names: &mut HashMap<u64, Vec<String>>,
    tip_names: &[String],
) -> TreeResult {
    let mut fingerprint_stack: Vec<u64> = Vec::new();
    let mut dist_stack: Vec<f64> = Vec::new();
    let mut tip_indices_stack: Vec<Vec<usize>> = Vec::new();
    let mut inner_nodes: Vec<InnerNodeInfo> = Vec::new();
    let mut split_fingerprints: Vec<u64> = Vec::new();
    let mut calibration_points: Vec<(f64, f64)> = Vec::new(); // (tip_date, dist_from_root)

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                let parent_dist = dist_stack.last().copied().unwrap_or(0.0);
                let dist = if dist_stack.is_empty() {
                    0.0 // root: ignore its branch_length
                } else {
                    parent_dist + node.branch_length
                };
                dist_stack.push(dist);

                if node.is_tip() {
                    let fp = *tip_fingerprints
                        .get(node.name.as_str())
                        .expect("Unknown tip name");
                    fingerprint_stack.push(fp);
                    let idx = tip_name_to_index[node.name.as_str()];
                    tip_indices_stack.push(vec![idx]);

                    if let Some(Some(date)) = tip_dates.get(node.name.as_str()) {
                        calibration_points.push((*date, dist));
                    }
                } else {
                    fingerprint_stack.push(0);
                    tip_indices_stack.push(Vec::new());
                }
            }
            TraversalAction::Exit => {
                let dist = dist_stack.pop().unwrap();
                let fp = fingerprint_stack.pop().unwrap();
                let indices = tip_indices_stack.pop().unwrap();

                if !node.is_tip() {
                    // Record clade tip names if this is a new fingerprint
                    if !clade_tip_names.contains_key(&fp) {
                        let mut names: Vec<String> =
                            indices.iter().map(|&i| tip_names[i].clone()).collect();
                        names.sort();
                        clade_tip_names.insert(fp, names);
                    }

                    inner_nodes.push(InnerNodeInfo {
                        fingerprint: fp,
                        dist_from_root: dist,
                    });
                    split_fingerprints.push(fp);
                }

                // Merge into parent
                if let Some(parent_fp) = fingerprint_stack.last_mut() {
                    *parent_fp ^= fp;
                }
                if let Some(parent_indices) = tip_indices_stack.last_mut() {
                    parent_indices.extend(indices);
                }
            }
        }
    }

    split_fingerprints.sort();

    let root_date = if calibration_points.is_empty() {
        None
    } else {
        let sum: f64 = calibration_points
            .iter()
            .map(|(tip_date, dist)| tip_date - dist)
            .sum();
        Some(sum / calibration_points.len() as f64)
    };

    TreeResult {
        sorted_split_fingerprints: split_fingerprints,
        inner_nodes,
        root_date,
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

fn compute_ess(sorted_splits: &[Vec<u64>]) -> f64 {
    let n = sorted_splits.len();
    if n <= 1 {
        return n as f64;
    }

    // All-pairs RF distance matrix
    let mut rf_dist = Array2::<f64>::zeros((n, n));
    for i in 0..n {
        eprintln!(
            "  Computing RF distances from tree {} of {}...",
            i + 1,
            n
        );
        for j in (i + 1)..n {
            let d = calc_rf_dist(&sorted_splits[i], &sorted_splits[j]) as f64;
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
    #[serde(skip_serializing_if = "Option::is_none")]
    mean_date: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    date_se: Option<f64>,
}

#[derive(Serialize)]
struct CladeResult {
    tips: Vec<String>,
    clade_size: usize,
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
    accumulators: HashMap<u64, CladeAccumulator>,
    ess: f64,
    num_total: usize,
    num_burnin: usize,
    num_post_burnin: usize,
    has_dates: bool,
}

fn process_file(
    file_path: &str,
    burnin: &BurninSpec,
    pool: &dyn Pool<NewickNode>,
    tip_fingerprints: &HashMap<String, u64>,
    tip_dates: &HashMap<String, Option<f64>>,
    tip_name_to_index: &HashMap<String, usize>,
    clade_tip_names: &mut HashMap<u64, Vec<String>>,
    tip_names: &[String],
) -> Result<FileProcessingResult, Box<dyn Error>> {
    eprintln!("Parsing {}...", file_path);
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let trees = NexusReader::parse(reader, 0, 1, pool)?;

    let num_total = trees.len();
    assert!(num_total >= 1, "No trees in {}", file_path);

    let num_burnin = burnin.first_sample_idx(num_total);
    let num_post_burnin = num_total - num_burnin;

    eprintln!(
        "  {} trees total, {} burn-in, {} post-burn-in",
        num_total, num_burnin, num_post_burnin
    );

    let mut accumulators: HashMap<u64, CladeAccumulator> = HashMap::new();
    let mut all_sorted_splits: Vec<Vec<u64>> = Vec::new();
    let mut has_dates = false;

    for (i, (_state, tree)) in trees.into_iter().enumerate() {
        if i < num_burnin {
            continue;
        }

        let result = process_tree(
            &tree,
            tip_fingerprints,
            tip_dates,
            tip_name_to_index,
            clade_tip_names,
            tip_names,
        );

        if result.root_date.is_some() {
            has_dates = true;
        }

        // Accumulate per-clade statistics
        for info in &result.inner_nodes {
            let acc = accumulators.entry(info.fingerprint).or_default();
            acc.count += 1;

            if let Some(root_date) = result.root_date {
                let node_date = root_date + info.dist_from_root;
                acc.date_sum += node_date;
                acc.date_sq_sum += node_date * node_date;
            }
        }

        all_sorted_splits.push(result.sorted_split_fingerprints);
    }

    eprintln!("Computing tree ESS for {}...", file_path);
    let ess = compute_ess(&all_sorted_splits);
    eprintln!("  ESS = {:.1}", ess);

    Ok(FileProcessingResult {
        accumulators,
        ess,
        num_total,
        num_burnin,
        num_post_burnin,
        has_dates,
    })
}

fn make_clade_file_stats(
    acc: Option<&CladeAccumulator>,
    num_post_burnin: usize,
    ess: f64,
    has_dates: bool,
) -> CladeFileStats {
    match acc {
        None | Some(CladeAccumulator { count: 0, .. }) => CladeFileStats {
            support: 0.0,
            support_se: 0.0,
            mean_date: if has_dates { Some(f64::NAN) } else { None },
            date_se: if has_dates { Some(f64::NAN) } else { None },
        },
        Some(acc) => {
            let p = acc.count as f64 / num_post_burnin as f64;
            let support_se = (p * (1.0 - p) / ess).sqrt();

            let (mean_date, date_se) = if has_dates && acc.count > 0 {
                let mean = acc.date_sum / acc.count as f64;
                let variance = (acc.date_sq_sum / acc.count as f64) - mean * mean;
                let sd = variance.max(0.0).sqrt();
                let ess_for_height = ess.min(acc.count as f64);
                let se = sd / ess_for_height.sqrt();
                (Some(mean), Some(se))
            } else {
                (None, None)
            };

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

    let burnin_spec = if let Some(pct) = args.burnin_pct {
        assert!((0.0..=100.0).contains(&pct));
        BurninSpec::Fract(pct / 100.0)
    } else if let Some(burnin_trees) = args.burnin_trees {
        BurninSpec::Trees(burnin_trees)
    } else {
        BurninSpec::Fract(0.10)
    };

    // Parse the first tree from file_a to get tip names and assign fingerprints
    eprintln!("Reading tip names from {}...", args.file_a);
    let file = File::open(&args.file_a)?;
    let reader = BufReader::new(file);
    let first_trees = NexusReader::parse(reader, 0, 1, &pool)?;
    assert!(!first_trees.is_empty(), "No trees in {}", args.file_a);

    let tip_fingerprints = assign_tip_fingerprints(&first_trees[0].1);

    // Build tip name list and index map
    let mut tip_names: Vec<String> = Vec::new();
    let mut tip_name_to_index: HashMap<String, usize> = HashMap::new();
    for node_ref in first_trees[0].1.any_order_iter() {
        let node = node_ref.borrow();
        if node.is_tip() {
            let idx = tip_names.len();
            tip_name_to_index.insert(node.name.clone(), idx);
            tip_names.push(node.name.clone());
        }
    }
    tip_names.sort();
    // Rebuild index map after sorting
    tip_name_to_index.clear();
    for (i, name) in tip_names.iter().enumerate() {
        tip_name_to_index.insert(name.clone(), i);
    }

    let num_tips = tip_names.len();
    eprintln!("  {} tips found", num_tips);

    // Build tip date map
    let mut tip_dates: HashMap<String, Option<f64>> = HashMap::new();
    let mut num_dated = 0;
    for name in &tip_names {
        let date = parse_tip_date(name);
        if date.is_some() {
            num_dated += 1;
        }
        tip_dates.insert(name.clone(), date);
    }
    eprintln!("  {} tips with full dates", num_dated);
    if num_dated == 0 {
        eprintln!("WARNING: No tips with full YYYY-MM-DD dates found. Date fields will be omitted.");
    }

    drop(first_trees); // Free memory before processing

    // Shared clade tip name map (fingerprint -> sorted tip names)
    let mut clade_tip_names: HashMap<u64, Vec<String>> = HashMap::new();

    // Process both files
    let result_a = process_file(
        &args.file_a,
        &burnin_spec,
        &pool,
        &tip_fingerprints,
        &tip_dates,
        &tip_name_to_index,
        &mut clade_tip_names,
        &tip_names,
    )?;
    let result_b = process_file(
        &args.file_b,
        &burnin_spec,
        &pool,
        &tip_fingerprints,
        &tip_dates,
        &tip_name_to_index,
        &mut clade_tip_names,
        &tip_names,
    )?;

    let has_dates = result_a.has_dates || result_b.has_dates;

    // Collect all clade fingerprints
    let mut all_fingerprints: Vec<u64> = clade_tip_names.keys().copied().collect();
    all_fingerprints.sort();

    // Build output, applying min-support filter
    eprintln!(
        "Building output ({} clades total, min-support={})...",
        all_fingerprints.len(),
        args.min_support
    );
    let mut clades: Vec<CladeResult> = Vec::new();
    for fp in &all_fingerprints {
        let stats_a = make_clade_file_stats(
            result_a.accumulators.get(fp),
            result_a.num_post_burnin,
            result_a.ess,
            has_dates,
        );
        let stats_b = make_clade_file_stats(
            result_b.accumulators.get(fp),
            result_b.num_post_burnin,
            result_b.ess,
            has_dates,
        );

        if stats_a.support < args.min_support && stats_b.support < args.min_support {
            continue;
        }

        let tips = clade_tip_names.get(fp).unwrap().clone();
        let clade_size = tips.len();

        clades.push(CladeResult {
            tips,
            clade_size,
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

        let pool = TestPool::new();
        let tree = NewickTree::new(pool.alloc(NewickNode::inner_node(
            "root",
            0.0,
            vec![
                pool.alloc(NewickNode::inner_node(
                    "",
                    1.0,
                    vec![
                        pool.alloc(NewickNode::leaf("A|2020-01-01", 0.5)),
                        pool.alloc(NewickNode::leaf("B|2020-07-01", 0.5)),
                    ],
                )),
                pool.alloc(NewickNode::leaf("C|2021-01-01", 2.0)),
            ],
        )));

        let tip_names = vec![
            "A|2020-01-01".to_string(),
            "B|2020-07-01".to_string(),
            "C|2021-01-01".to_string(),
        ];
        let tip_name_to_index: HashMap<String, usize> = tip_names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        let tip_fingerprints: HashMap<String, u64> = HashMap::from([
            ("A|2020-01-01".to_string(), 0b001),
            ("B|2020-07-01".to_string(), 0b010),
            ("C|2021-01-01".to_string(), 0b100),
        ]);
        let tip_dates: HashMap<String, Option<f64>> = tip_names
            .iter()
            .map(|n| (n.clone(), parse_tip_date(n)))
            .collect();

        let mut clade_tip_names: HashMap<u64, Vec<String>> = HashMap::new();

        let result = process_tree(
            &tree,
            &tip_fingerprints,
            &tip_dates,
            &tip_name_to_index,
            &mut clade_tip_names,
            &tip_names,
        );

        assert!(result.root_date.is_some());
        // 2 inner nodes: {A,B} and the root {A,B,C}
        assert_eq!(result.inner_nodes.len(), 2);
        assert_eq!(result.sorted_split_fingerprints.len(), 2);
        // Should have discovered clades
        assert!(!clade_tip_names.is_empty());
    }
}
