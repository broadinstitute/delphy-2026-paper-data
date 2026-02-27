use ndarray::{Array1, ArrayView2};

pub struct FrechetEssResult {
    pub rho_s: Vec<f64>,
    pub auto_correlation_time: f64,
    pub frechet_ess: f64,
}

/// Calculate the Frechet correlation ESS from an all-pairs distance matrix.
///
/// Implements the method from Magee et al, "How Trustworthy Is Your Tree?
/// Bayesian Phylogenetic Effective Sample Size Through the Lens of Monte
/// Carlo Error", Bayesian Analysis 19(2), 565--593 (2024).
///
/// We simplify their expressions because Var(tau_t) and Var(tau_{t+s}) are
/// statistically identical (why should the variance at different points in
/// the chain be different when averaged over all possible chains?), and we
/// can thus use all pairwise distances to get a slightly better estimate
/// for them (Eq (15) instead of the equations following Eq (18)).
/// With that simplification, we get,
///
///   rho_s = 1 - (1/2) E[Delta_s^2] / Var[tau],
///
/// where
///
///   Var[tau] = 1/(n (n-1)) sum_{i=1}^n sum_{j=i+1}^n d^2(x_i, x_j)
///   E[Delta_s^2] = 1/(n-s) sum_{i=1}^{n-s} d^2(x_i, x_{i+s})
pub fn calc_frechet_ess(d_ij: &ArrayView2<f64>) -> FrechetEssResult {
    let n = d_ij.nrows();

    if n <= 1 {
        return FrechetEssResult {
            rho_s: vec![],
            auto_correlation_time: 1.0,
            frechet_ess: n as f64,
        };
    }

    // Var[tau] = 1/(n(n-1)) sum_{i<j} d^2(x_i, x_j)
    let mut var_tau = 0.0;
    for i in 0..n {
        for j in (i + 1)..n {
            let d = d_ij[(i, j)];
            var_tau += d * d;
        }
    }
    var_tau /= (n as f64) * ((n - 1) as f64);

    if var_tau == 0.0 {
        return FrechetEssResult {
            rho_s: vec![],
            auto_correlation_time: 1.0,
            frechet_ess: n as f64,
        };
    }

    // E[Delta_s^2] at each lag s
    let mut mean_delta_s_2 = Array1::<f64>::zeros(n);
    for s in 0..n {
        for i in 0..(n - s) {
            let d = d_ij[(i, i + s)];
            mean_delta_s_2[s] += d * d;
        }
        mean_delta_s_2[s] /= (n - s) as f64;
    }

    let rho_s = 1.0 - 0.5 * mean_delta_s_2 / var_tau;

    // Calculate auto-correlation time as 1 + 2 * sum_{s=1}^{smax-1} rho_s[s],
    // where smax is the first s for which rho_s[s-1] + rho_s[s] <= 0
    // (i.e., rho_s[s] drops below 0, and it's not just a fluke)
    let mut smax = 1;
    while (smax < n) && ((rho_s[smax - 1] + rho_s[smax]) > 0.0) {
        smax += 1;
    }
    let mut act = 1.0;
    for s in 1..smax {
        act += 2.0 * rho_s[s];
    }
    let ess = n as f64 / act;

    FrechetEssResult {
        rho_s: rho_s.to_vec(),
        auto_correlation_time: act,
        frechet_ess: ess,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::array;

    #[test]
    fn simple_3x3() {
        // d = [[0, 2, 4],
        //      [2, 0, 2],
        //      [4, 2, 0]]
        //
        // var_tau = (4+16+4) / (3*2) = 24/6 = 4
        // mean_delta_s_2[0] = (0+0+0)/3 = 0
        // mean_delta_s_2[1] = (4+4)/2 = 4
        // mean_delta_s_2[2] = 16/1 = 16
        // rho_s = [1.0, 0.5, -1.0]
        // smax = 2 (rho_s[1]+rho_s[2] = -0.5 <= 0)
        // ACT = 1 + 2*0.5 = 2.0
        // ESS = 3/2 = 1.5
        let d = array![[0.0, 2.0, 4.0],
                        [2.0, 0.0, 2.0],
                        [4.0, 2.0, 0.0]];

        let result = calc_frechet_ess(&d.view());

        assert_eq!(result.rho_s, vec![1.0, 0.5, -1.0]);
        assert_eq!(result.auto_correlation_time, 2.0);
        assert_eq!(result.frechet_ess, 1.5);
    }
}
