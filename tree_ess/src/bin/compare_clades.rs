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
use ndarray::Array2;
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use serde::Serialize;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io;
use tree_ess::burnin::BurninSpec;
use tree_ess::clades::{calc_clade_times_from_root, catalog_tree_clades, assign_tip_fps, calc_rf_dist, tips_in_clade, CladeDefinition, CladeFp, CladeMap};
use tree_ess::dates::find_root_date;
use tree_ess::ess::calc_frechet_ess;
use tree_ess::newick::NewickTree;
use tree_ess::nexus_reader::NexusReader;
use tree_ess::refs::AllocPool;
use tree_ess::trees::{NodeLike, TreeLike};

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

// -- Per-tree processing --

struct TreeResult {
    sorted_clade_fps: Vec<CladeFp>,
    clade_fps_2_dates: HashMap<CladeFp, f64>,
}

fn process_tree(
    tree: &NewickTree,
    tip_fps: &HashMap<String, CladeFp>,
    clade_map: &mut CladeMap,
) -> TreeResult {

    let sorted_clade_fps = catalog_tree_clades(tree, tip_fps, clade_map, true);
    let clade_times_from_root = calc_clade_times_from_root(tree, tip_fps);

    let root_date = find_root_date(tree)
        .expect("No tip found with an exact date!  Cannot deduce root date!");

    let clade_fps_2_dates: HashMap<CladeFp, f64> = clade_times_from_root
        .into_iter()
        .map(|(fp, time_from_root)| (fp, root_date + time_from_root))
        .collect();

    TreeResult {
        sorted_clade_fps,
        clade_fps_2_dates,
    }
}

// -- Clade accumulator --

#[derive(Default)]
struct CladeAccumulator {
    count: usize,
    date_sum: f64,
    date_sq_sum: f64,
}

// -- Tree ESS (Frechet correlation ESS over RF distance matrix) --

fn calc_tree_ess(sorted_clade_fps: &[Vec<CladeFp>]) -> f64 {
    let n = sorted_clade_fps.len();

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

    calc_frechet_ess(&rf_dist.view()).frechet_ess
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

        let result = process_tree(&tree, tip_fps, clade_map);

        // Accumulate per-clade statistics
        for (&clade_fp, &date) in &result.clade_fps_2_dates {
            let acc = accumulators.entry(clade_fp).or_default();
            acc.count += 1;
            acc.date_sum += date;
            acc.date_sq_sum += date * date;
        }

        all_sorted_clade_fps.push(result.sorted_clade_fps);
    }

    eprintln!("Computing tree ESS for {}...", file_path);
    let ess = calc_tree_ess(&all_sorted_clade_fps);
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

    // Clade definitions will accumulate here
    let mut clade_map: CladeMap = HashMap::new();

    // Process both files
    let result_a = process_trees(
        &args.file_a,
        &file_a_trees,
        &burnin_spec,
        &tip_fps,
        &mut clade_map,
    )?;
    let result_b = process_trees(
        &args.file_b,
        &file_b_trees,
        &burnin_spec,
        &tip_fps,
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
    use tree_ess::clades::CladeDefinition::{InnerNodeClade, TipClade};
    use itertools::Itertools;
    use tree_ess::newick::NewickNode;
    use tree_ess::refs::Pool;

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
        let mut clade_map: CladeMap = HashMap::new();

        let result = process_tree(&tree, &tip_fps, &mut clade_map);

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
                    }
                ),
                (
                    fp_abc,
                    InnerNodeClade {
                        subclade1: fp_ab,
                        subclade2: fp_c,
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
