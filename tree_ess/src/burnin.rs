use std::cmp;

#[derive(Debug)]
pub enum BurninSpec {
    Fract(f64),
    Trees(usize),
}

impl BurninSpec {
    /// Construct from CLI options.  `pct` and `trees` are mutually
    /// exclusive (enforced by clap's `group`).  If neither is given,
    /// defaults to 10% burn-in.
    pub fn from_options(pct: Option<f64>, trees: Option<usize>) -> Self {
        if let Some(pct) = pct {
            assert!((0.0..=100.0).contains(&pct));
            BurninSpec::Fract(pct / 100.0)
        } else if let Some(trees) = trees {
            BurninSpec::Trees(trees)
        } else {
            BurninSpec::Fract(0.10)
        }
    }

    pub fn first_sample_idx(&self, num_trees: usize) -> usize {
        match *self {
            BurninSpec::Fract(pct) => {
                assert!((0.0..=1.0).contains(&pct));
                (num_trees as f64 * pct).floor() as usize
            }
            BurninSpec::Trees(burnin_trees) => cmp::min(num_trees, burnin_trees),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn from_options_default_is_10_percent() {
        let spec = BurninSpec::from_options(None, None);
        assert_eq!(spec.first_sample_idx(100), 10);
        assert_eq!(spec.first_sample_idx(1000), 100);
    }

    #[test]
    fn from_options_pct() {
        let spec = BurninSpec::from_options(Some(25.0), None);
        assert_eq!(spec.first_sample_idx(100), 25);
        assert_eq!(spec.first_sample_idx(200), 50);
    }

    #[test]
    fn from_options_trees() {
        let spec = BurninSpec::from_options(None, Some(50));
        assert_eq!(spec.first_sample_idx(100), 50);
        assert_eq!(spec.first_sample_idx(200), 50);
    }

    #[test]
    fn trees_capped_at_num_trees() {
        let spec = BurninSpec::from_options(None, Some(500));
        assert_eq!(spec.first_sample_idx(100), 100);
    }

    #[test]
    fn zero_burnin_pct() {
        let spec = BurninSpec::from_options(Some(0.0), None);
        assert_eq!(spec.first_sample_idx(100), 0);
    }

    #[test]
    fn full_burnin_pct() {
        let spec = BurninSpec::from_options(Some(100.0), None);
        assert_eq!(spec.first_sample_idx(100), 100);
    }

    #[test]
    fn fractional_burnin_floors() {
        // 10% of 15 = 1.5, should floor to 1
        let spec = BurninSpec::from_options(Some(10.0), None);
        assert_eq!(spec.first_sample_idx(15), 1);
    }

    #[test]
    #[should_panic]
    fn pct_over_100_panics() {
        BurninSpec::from_options(Some(101.0), None);
    }

    #[test]
    #[should_panic]
    fn negative_pct_panics() {
        BurninSpec::from_options(Some(-1.0), None);
    }
}
