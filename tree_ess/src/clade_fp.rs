use rand::Rng;
use serde::Serialize;

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
}
