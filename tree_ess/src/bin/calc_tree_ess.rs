//! Calculate a "tree ESS" for a BEAST-like .trees file using the frechetCorrelationESS
//! metric from Magee, et al, "How Trustworthy Is Your Tree? Bayesian Phylogenetic Effective
//! Sample Size Through the Lens of Monte Carlo Error", Bayesian Analysis 19(2), 565--593 (2024).

extern crate tree_ess;

use clap::Parser;
use ndarray::{s, Array2, ArrayView2};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use serde::Serialize;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::io;
use tree_ess::burnin::BurninSpec;
use tree_ess::clades::{assign_tip_fps, calc_rf_dist, CladeFp};
use tree_ess::ess::calc_frechet_ess;
use tree_ess::newick::NewickTree;
use tree_ess::nexus_reader::NexusReader;
use tree_ess::refs::AllocPool;
use tree_ess::trees::NodeLike;
use tree_ess::trees::{TraversalAction, TreeLike};

#[derive(Parser)]
#[command(version, about, long_about = None)]
#[command(name = "calc_tree_ess")]
struct Args {
    /// BEAST-like .trees files over which to calculate tree ESS
    #[arg(required = true)]
    files: Vec<String>,

    /// Percentage of initial trees to discard as burn-in
    #[arg(long, group = "burn-in")]
    burnin_pct: Option<f64>,

    /// Number of initial trees to discard as burn-in
    #[arg(long, group = "burn-in")]
    burnin_trees: Option<usize>,

    /// Enable compact output format (don't output whole Robison-Foulds distance matrix)
    #[arg(short, long)]
    compact: bool,
}

struct Chain {
    file_path: String,
    samples: Vec<Sample>,
    results: Option<ChainResults>,
}
struct Sample {
    global_sample_num: usize,
    state: u64,
    //tree: NewickTree,
    sorted_split_fps: Vec<CladeFp>,
}

#[derive(Serialize)]
struct FinalResult {
    chain_results: Vec<ChainResults>,
    all_d_ij: Option<Vec<Vec<f64>>>,
}

#[derive(Serialize)]
struct ChainResults {
    num_samples: u64,
    burn_in_samples: usize,
    d_ij: Option<Vec<Vec<f64>>>,
    rho_s: Vec<f64>,
    auto_correlation_time: f64,
    effective_sample_size: f64,
}

fn array2_to_vec_vec<T: Clone>(arr: &ArrayView2<T>) -> Vec<Vec<T>> {
    arr.rows().into_iter().map(|row| row.to_vec()).collect()
}

fn main() -> Result<(), Box<dyn Error>> {
    let pool = AllocPool::new();

    let args = Args::parse();
    let burnin_spec = BurninSpec::from_options(args.burnin_pct, args.burnin_trees);

    let mut rng = Pcg64Mcg::seed_from_u64(1);

    let mut chains = vec![];
    let mut next_global_sample_num = 0;
    let mut tip_name_2_fp_opt = None;
    for file_path in args.files {
        let file = File::open(&file_path)?;
        let reader = BufReader::new(file);
        let trees = NexusReader::parse(reader, 0, 1, &pool)?;

        assert!(
            trees.len() >= 1,
            "Can't calculate an RS distance matrix for {} trees in {}",
            trees.len(),
            file_path
        );

        // Assign fingerprints to tips in first tree (the same tips should appear in every tree!)
        let tip_name_2_fp =
            tip_name_2_fp_opt.get_or_insert_with(|| assign_tip_fps(&trees[0].1, &mut rng));

        let mut samples = vec![];
        for (state, tree) in trees.into_iter() {
            // Now list out all (rooted) splits in each tree.
            // Each branch of the tree splits the nodes of the tree into sets {X,Y}, and the root node is
            // in exactly one of those, wlog Y.  We represent the split by the fingerprint of the
            // nodes in X (all inner nodes are unlabeled).  Since the splits {{<single tip>}, {<rest>}}
            // are in every tree, we don't materialize those
            let sorted_splits = calc_sorted_splits(&tree, &tip_name_2_fp);

            samples.push(Sample {
                global_sample_num: next_global_sample_num,
                state,
                //tree,
                sorted_split_fps: sorted_splits,
            });

            next_global_sample_num += 1;
        }

        chains.push(Chain {
            file_path,
            samples,
            results: None,
        });
    }
    let total_samples = next_global_sample_num;

    assert!(
        chains.len() >= 1,
        "Can't calculate an RS distance matrix for 0 chains"
    );

    // Calculate all-to-all RF distance matrix
    let mut all_rf_dist_matrix = Array2::<f64>::zeros((total_samples, total_samples));

    for (a, chain_a) in chains.iter().enumerate() {
        let n_a = chain_a.samples.len();
        for (b, chain_b) in chains.iter().enumerate() {
            let n_b = chain_b.samples.len();

            if b < a {
                continue;
            }

            // And here we go!  Use f64 so that later calculations are easier to write out
            for i in 0..n_a {
                let sample_i = &chain_a.samples[i];
                eprintln!(
                    "Calculating distances from state {} state {} to all states in {}",
                    chain_a.file_path, sample_i.state, chain_b.file_path
                );

                let ii = sample_i.global_sample_num;
                let sorted_split_fps_i = &sample_i.sorted_split_fps;
                for j in (i + 1)..n_b {
                    let sample_j = &chain_b.samples[j];
                    let jj = sample_j.global_sample_num;
                    let sorted_split_fps_j = &sample_j.sorted_split_fps;

                    let d_u64 =
                        calc_rf_dist(sorted_split_fps_i, sorted_split_fps_j);
                    let d = d_u64 as f64;

                    all_rf_dist_matrix[(ii, jj)] = d;
                    all_rf_dist_matrix[(jj, ii)] = d;
                }
            }
        }
    }

    for chain in chains.iter_mut() {
        let n = chain.samples.len();

        // Extract the intra-chain RF distance matrices
        let mut rf_dist_matrix = Array2::<f64>::zeros((n, n));
        for (i, sample_i) in chain.samples.iter().enumerate() {
            let ii = sample_i.global_sample_num;
            for (j, sample_j) in chain.samples.iter().enumerate() {
                let jj = sample_j.global_sample_num;

                rf_dist_matrix[(i, j)] = all_rf_dist_matrix[(ii, jj)];
            }
        }
        let rf_dist_matrix = rf_dist_matrix; // End of mutation

        // Burn-in
        let first_prod_sample_idx = burnin_spec.first_sample_idx(n);
        let prod_d_ij =
            rf_dist_matrix.slice(s![first_prod_sample_idx..n, first_prod_sample_idx..n]);
        let prod_n = n - first_prod_sample_idx;

        let frechet_ess_result = calc_frechet_ess(&prod_d_ij);

        chain.results = Some(ChainResults {
            num_samples: chain.samples.len() as u64,
            burn_in_samples: first_prod_sample_idx,
            d_ij: Some(array2_to_vec_vec(&rf_dist_matrix.view())),
            rho_s: frechet_ess_result.rho_s,
            auto_correlation_time: frechet_ess_result.auto_correlation_time,
            effective_sample_size: frechet_ess_result.frechet_ess,
        });
        eprintln!(
            "For {}: N = {}, ACT = {}, ESS = {}",
            chain.file_path, prod_n, frechet_ess_result.auto_correlation_time, frechet_ess_result.frechet_ess
        );
    }

    // Pedestrian json output
    eprintln!("Writing output...");
    let final_result = FinalResult {

        chain_results: chains
            .into_iter()
            .map(|mut chain| {
                let results = chain.results.take().unwrap();
                if args.compact {
                    ChainResults {
                        d_ij: None,
                        ..results
                    }
                } else {
                    results
                }
            })
            .collect(),

        all_d_ij: if args.compact {
            None
        } else {
            Some(array2_to_vec_vec(&all_rf_dist_matrix.view()))
        },
    };
    serde_json::to_writer_pretty(io::stdout(), &final_result)?;

    eprintln!("Done!");

    Ok(())
}

fn calc_sorted_splits(
    tree: &NewickTree,
    tip_name_2_fp: &HashMap<String, CladeFp>,
) -> Vec<CladeFp> {
    let mut partial_fps: Vec<CladeFp> = Vec::new();
    let mut result: Vec<CladeFp> = Vec::new();
    for (action, node_ref) in tree.traversal_iter() {
        match action {
            TraversalAction::Enter => {
                if node_ref.borrow().is_tip() {
                    partial_fps.push(
                        *tip_name_2_fp
                            .get(node_ref.borrow().name.as_str())
                            .unwrap(),
                    );
                } else {
                    partial_fps.push(CladeFp::empty());
                }
            }
            TraversalAction::Exit => {
                let fp = partial_fps.pop().expect("Improper nesting?");
                if !node_ref.borrow().is_tip() {
                    result.push(fp);
                }
                if let Some(parent_fp) = partial_fps.last_mut() {
                    *parent_fp = parent_fp.union(&fp); // fp(node) = XOR_{children i} fp(i)
                }
            }
        }
    }
    result.sort();
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use tree_ess::newick::NewickNode;
    use tree_ess::refs::{Pool, TestPool};

    #[test]
    fn calc_sorted_splits_test() {
        let pool = TestPool::new();
        let tree = NewickTree::new(pool.alloc(NewickNode::inner_node(
            "root",
            0.0,
            vec![
                pool.alloc(NewickNode::inner_node(
                    "X",
                    0.0,
                    vec![
                        pool.alloc(NewickNode::leaf("A", 0.0)),
                        pool.alloc(NewickNode::leaf("B", 0.0)),
                    ],
                )),
                pool.alloc(NewickNode::inner_node(
                    "X",
                    0.0,
                    vec![
                        pool.alloc(NewickNode::leaf("C", 0.0)),
                        pool.alloc(NewickNode::leaf("D", 0.0)),
                    ],
                )),
            ],
        )));
        let tip_name_2_fp = HashMap::from([
            (String::from("A"), CladeFp::new(0b0001)),
            (String::from("B"), CladeFp::new(0b0010)),
            (String::from("C"), CladeFp::new(0b0100)),
            (String::from("D"), CladeFp::new(0b1000)),
        ]);

        let splits = calc_sorted_splits(&tree, &tip_name_2_fp);

        assert_eq!(
            splits,
            vec![
                CladeFp::new(0b0011), // {A,B} || {C,D,root}
                CladeFp::new(0b1100), // {C,D} || {A,B,root}
                CladeFp::new(0b1111), // {A,B,C,D} || {root}
            ]
        );
    }
}
