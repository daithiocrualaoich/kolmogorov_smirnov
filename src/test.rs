//! Two Sample Kolmogorov-Smirnov Test

use std::cmp::{min, Ord, Ordering};

/// Two sample test result.
pub struct TestResult {
    pub is_rejected: bool,
    pub statistic: f64,
    pub critical_value: f64,
    pub confidence: f64,
}

/// Perform a two sample Kolmogorov-Smirnov test on given samples.
///
/// The samples currently must have length > 12 elements for the test to be
/// valid. Also, only the 0.95 confidence level is supported initially.
///
/// # Panics
///
/// There are assertion panics if either sequence has <= 12 elements or if
/// confidence is not 0.95.
///
/// # Examples
///
/// ```
/// extern crate kolmogorov_smirnov as ks;
///
/// let xs = vec!(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
/// let ys = vec!(12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
/// let confidence = 0.95;
///
/// let result = ks::test(&xs, &ys, confidence);
///
/// if result.is_rejected {
///     println!("{:?} and {:?} are not from the same distribution with confidence {}.",
///       xs, ys, confidence);
/// }
/// ```
pub fn test<T: Ord + Clone>(xs: &[T], ys: &[T], confidence: f64) -> TestResult {
    assert!(xs.len() > 0 && ys.len() > 0);
    assert!(0.0 < confidence && confidence < 1.0);

    // Only support samples of size > 12 initially.
    assert!(xs.len() > 12 && ys.len() > 12);

    // Only support confidence == 0.95 initially.
    assert_eq!(confidence, 0.95);

    let statistic = calculate_statistic(xs, ys);
    let critical_value = calculate_critical_value(xs.len(), ys.len(), confidence);
    let is_rejected = statistic > critical_value;

    TestResult {
        is_rejected: is_rejected,
        statistic: statistic,
        critical_value: critical_value,
        confidence: confidence,
    }
}

/// Wrapper type for f64 to implement Ord and make usable with test.
#[derive(PartialEq, Clone)]
struct OrderableF64 {
    val: f64,
}

impl OrderableF64 {
    fn new(val: f64) -> OrderableF64 {
        OrderableF64 { val: val }
    }
}

impl Eq for OrderableF64 {}

impl PartialOrd for OrderableF64 {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        self.val.partial_cmp(&other.val)
    }
}

impl Ord for OrderableF64 {
    fn cmp(&self, other: &Self) -> Ordering {
        self.val.partial_cmp(&other.val).unwrap()
    }
}

/// Perform a two sample Kolmogorov-Smirnov test on given f64 samples.
///
/// This is necessary because f64 does not implement Ord in Rust as some
/// elements are incomparable, e.g. NaN. This function wraps the f64s in
/// implementation of Ord which panics on incomparable elements.
///
/// The samples currently must have length > 12 elements for the test to be
/// valid. Also, only the 0.95 confidence level is supported initially.
///
/// # Panics
///
/// There are assertion panics if either sequence has <= 12 elements or if
/// confidence is not 0.95.
///
/// If any of the f64 elements in the input samples are unorderable, e.g. NaN.
///
/// # Examples
///
/// ```
/// extern crate kolmogorov_smirnov as ks;
///
/// let xs = vec!(0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0);
/// let ys = vec!(12.0, 11.0, 10.0, 9.0, 8.0, 7.0, 6.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0);
/// let confidence = 0.95;
///
/// let result = ks::test_f64(&xs, &ys, confidence);
///
/// if result.is_rejected {
///     println!("{:?} and {:?} are not from the same distribution with confidence {}.",
///       xs, ys, confidence);
/// }
/// ```
pub fn test_f64(xs: &[f64], ys: &[f64], confidence: f64) -> TestResult {
    let xs: Vec<OrderableF64> = xs.iter().map(|&f| OrderableF64::new(f)).collect();
    let ys: Vec<OrderableF64> = ys.iter().map(|&f| OrderableF64::new(f)).collect();

    test(&xs, &ys, confidence)
}

/// Calculate the critical value for the two sample Kolmogorov-Smirnov test.
fn calculate_critical_value(n1: usize, n2: usize, confidence: f64) -> f64 {
    assert!(n1 > 0 && n2 > 0);
    assert!(0.0 < confidence && confidence < 1.0);

    // Only support samples of size > 12 initially.
    assert!(n1 > 12 && n2 > 12);

    // Only support confidence == 0.95 initially.
    assert_eq!(confidence, 0.95);

    let n1 = n1 as f64;
    let n2 = n2 as f64;

    let factor = (n1 + n2) / (n1 * n2);
    1.36 * factor.sqrt()
}

/// Calculate the test statistic for the two sample Kolmogorov-Smirnov test.
///
/// The test statistic is the maximum vertical distance between the ECDFs of
/// the two samples.
fn calculate_statistic<T: Ord + Clone>(xs: &[T], ys: &[T]) -> f64 {
    let n = xs.len();
    let m = ys.len();

    assert!(n > 0 && m > 0);

    let mut xs = xs.to_vec();
    let mut ys = ys.to_vec();

    // xs and ys must be sorted for the stepwise ECDF calculations to work.
    xs.sort();
    ys.sort();

    // The current value testing for ECDF difference. Sweeps up through elements
    // present in xs and ys.
    let mut current: &T;

    // i, j index the first values in xs and ys that are greater than current.
    let mut i = 0;
    let mut j = 0;

    // ecdf_xs, ecdf_ys always hold the ECDF(current) of xs and ys.
    let mut ecdf_xs = 0.0;
    let mut ecdf_ys = 0.0;

    // The test statistic value computed over values <= current.
    let mut statistic = 0.0;

    while i < n && j < m {
        // Advance i through duplicate samples in xs.
        let x_i = &xs[i];
        while i + 1 < n && *x_i == xs[i + 1] {
            i += 1;
        }

        // Advance j through duplicate samples in ys.
        let y_j = &ys[j];
        while j + 1 < m && *y_j == ys[j + 1] {
            j += 1;
        }

        // Step to the next sample value in the ECDF sweep from low to high.
        current = min(x_i, y_j);

        // Update invariant conditions for i, j, ecdf_xs, and ecdf_ys.
        if current == x_i {
            ecdf_xs = (i + 1) as f64 / n as f64;
            i += 1;
        }
        if current == y_j {
            ecdf_ys = (j + 1) as f64 / m as f64;
            j += 1;
        }

        // Update invariant conditions for the test statistic.
        let diff = (ecdf_xs - ecdf_ys).abs();
        if diff > statistic {
            statistic = diff;
        }
    }

    // Don't need to walk the rest of the samples because one of the ecdfs is
    // already one and the other will be increasing up to one. This means the
    // difference will be monotonically decreasing, so we have our test
    // statistic value already.

    statistic
}


#[cfg(test)]
mod tests {
    extern crate quickcheck;
    extern crate rand;

    use self::quickcheck::{Arbitrary, Gen, QuickCheck, Testable, StdGen};
    use self::rand::Rng;
    use std::cmp;
    use std::usize;

    use super::test;
    use ecdf::Ecdf;

    const EPSILON: f64 = 1e-10;

    fn check<A: Testable>(f: A) {
        // Need - 1 to ensure space for creating non-overlapping samples.
        let g = StdGen::new(rand::thread_rng(), usize::MAX - 1);
        QuickCheck::new().gen(g).quickcheck(f);
    }

    /// Wrapper for generating sample data with QuickCheck.
    ///
    /// Samples must be sequences of u64 values with more than 12 elements.
    #[derive(Debug, Clone)]
    struct Samples {
        vec: Vec<u64>,
    }

    impl Samples {
        fn min(&self) -> u64 {
            let &min = self.vec.iter().min().unwrap();
            min
        }

        fn max(&self) -> u64 {
            let &max = self.vec.iter().max().unwrap();
            max
        }

        fn shuffle(&mut self) {
            let mut rng = rand::thread_rng();
            rng.shuffle(&mut self.vec);
        }
    }

    impl Arbitrary for Samples {
        fn arbitrary<G: Gen>(g: &mut G) -> Samples {
            // Limit size of generated sample set to 1024
            let max = cmp::min(g.size(), 1024);

            let size = g.gen_range(13, max);
            let vec = (0..size).map(|_| u64::arbitrary(g)).collect();

            Samples { vec: vec }
        }

        fn shrink(&self) -> Box<Iterator<Item = Samples>> {
            let vec: Vec<u64> = self.vec.clone();
            let shrunk: Box<Iterator<Item = Vec<u64>>> = vec.shrink();

            Box::new(shrunk.filter(|v| v.len() > 12).map(|v| Samples { vec: v }))
        }
    }

    #[test]
    #[should_panic(expected="assertion failed: xs.len() > 0 && ys.len() > 0")]
    fn test_panics_on_empty_samples_set() {
        let xs: Vec<u64> = vec![];
        let ys: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        test(&xs, &ys, 0.95);
    }

    #[test]
    #[should_panic(expected="assertion failed: xs.len() > 0 && ys.len() > 0")]
    fn test_panics_on_empty_other_samples_set() {
        let xs: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        let ys: Vec<u64> = vec![];
        test(&xs, &ys, 0.95);
    }

    #[test]
    #[should_panic(expected="assertion failed: 0.0 < confidence && confidence < 1.0")]
    fn test_panics_on_confidence_leq_zero() {
        let xs: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        let ys: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        test(&xs, &ys, 0.0);
    }

    #[test]
    #[should_panic(expected="assertion failed: 0.0 < confidence && confidence < 1.0")]
    fn test_panics_on_confidence_geq_one() {
        let xs: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        let ys: Vec<u64> = vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12];
        test(&xs, &ys, 1.0);
    }

    #[test]
    fn test_is_rejected_if_test_statistic_greater_than_critical_value() {
        fn prop(xs: Samples, ys: Samples) -> bool {
            let result = test(&xs.vec, &ys.vec, 0.95);

            if result.is_rejected {
                result.statistic > result.critical_value
            } else {
                result.statistic <= result.critical_value
            }
        }

        check(prop as fn(Samples, Samples) -> bool);
    }

    /// Alternative calculation for the test statistic for the two sample
    /// Kolmogorov-Smirnov test. This simple implementation is used as a
    /// verification check against actual calculation used.
    fn calculate_statistic_alt<T: Ord + Clone>(xs: &[T], ys: &[T]) -> f64 {
        assert!(xs.len() > 0 && ys.len() > 0);

        let ecdf_xs = Ecdf::new(xs);
        let ecdf_ys = Ecdf::new(ys);

        let mut statistic = 0.0;

        for x in xs.iter() {
            let diff = (ecdf_xs.value(x.clone()) - ecdf_ys.value(x.clone())).abs();
            if diff > statistic {
                statistic = diff;
            }
        }

        for y in ys.iter() {
            let diff = (ecdf_xs.value(y.clone()) - ecdf_ys.value(y.clone())).abs();
            if diff > statistic {
                statistic = diff;
            }
        }

        statistic
    }

    #[test]
    fn test_calculate_statistic() {
        fn prop(xs: Samples, ys: Samples) -> bool {
            let result = test(&xs.vec, &ys.vec, 0.95);
            let actual = result.statistic;
            let expected = calculate_statistic_alt(&xs.vec, &ys.vec);

            actual == expected
        }

        check(prop as fn(Samples, Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_between_zero_and_one() {
        fn prop(xs: Samples, ys: Samples) -> bool {
            let result = test(&xs.vec, &ys.vec, 0.95);
            let actual = result.statistic;

            0.0 <= actual && actual <= 1.0
        }

        check(prop as fn(Samples, Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_zero_for_identical_samples() {
        fn prop(xs: Samples) -> bool {
            let ys = xs.clone();

            let result = test(&xs.vec, &ys.vec, 0.95);

            result.statistic == 0.0
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_zero_for_permuted_sample() {
        fn prop(xs: Samples) -> bool {
            let mut ys = xs.clone();
            ys.shuffle();

            let result = test(&xs.vec, &ys.vec, 0.95);

            result.statistic == 0.0
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_one_for_samples_with_no_overlap_in_support() {
        fn prop(xs: Samples) -> bool {
            let mut ys = xs.clone();

            // Shift ys so that ys.min > xs.max.
            let ys_min = xs.max() + 1;
            ys.vec = ys.vec.iter().map(|&y| cmp::max(y, ys_min)).collect();

            let result = test(&xs.vec, &ys.vec, 0.95);

            result.statistic == 1.0
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_one_half_for_sample_with_non_overlapping_in_support_replicate_added() {
        fn prop(xs: Samples) -> bool {
            let mut ys = xs.clone();

            // Shift ys so that ys.min > xs.max.
            let ys_min = xs.max() + 1;
            ys.vec = ys.vec.iter().map(|&y| cmp::max(y, ys_min)).collect();

            // Add all the original items back too.
            for &x in xs.vec.iter() {
                ys.vec.push(x);
            }

            let result = test(&xs.vec, &ys.vec, 0.95);

            result.statistic == 0.5
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_one_div_length_for_sample_with_additional_low_value() {
        fn prop(xs: Samples) -> bool {
            // Add a extra sample of early weight to ys.
            let min = xs.min();
            let mut ys = xs.clone();
            ys.vec.push(min - 1);

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = 1.0 / ys.vec.len() as f64;

            result.statistic == expected
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_one_div_length_for_sample_with_additional_high_value() {
        fn prop(xs: Samples) -> bool {
            // Add a extra sample of late weight to ys.
            let max = xs.max();
            let mut ys = xs.clone();
            ys.vec.push(max + 1);

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = 1.0 / ys.vec.len() as f64;

            (result.statistic - expected).abs() < EPSILON
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_one_div_length_for_sample_with_additional_low_and_high_values() {
        fn prop(xs: Samples) -> bool {
            // Add a extra sample of late weight to ys.
            let min = xs.min();
            let max = xs.max();

            let mut ys = xs.clone();

            ys.vec.push(min - 1);
            ys.vec.push(max + 1);

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = 1.0 / ys.vec.len() as f64;

            (result.statistic - expected).abs() < EPSILON
        }

        check(prop as fn(Samples) -> bool);
    }

    #[test]
    fn test_statistic_is_n_div_length_for_sample_with_additional_n_low_values() {
        fn prop(xs: Samples, n: u8) -> bool {
            // Add extra sample of early weight to ys.
            let min = xs.min();
            let mut ys = xs.clone();
            for j in 0..n {
                ys.vec.push(min - (j as u64) - 1);
            }

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = n as f64 / ys.vec.len() as f64;

            result.statistic == expected
        }

        check(prop as fn(Samples, u8) -> bool);
    }

    #[test]
    fn test_statistic_is_n_div_length_for_sample_with_additional_n_high_values() {
        fn prop(xs: Samples, n: u8) -> bool {
            // Add extra sample of early weight to ys.
            let max = xs.max();
            let mut ys = xs.clone();
            for j in 0..n {
                ys.vec.push(max + (j as u64) + 1);
            }

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = n as f64 / ys.vec.len() as f64;

            (result.statistic - expected).abs() < EPSILON
        }

        check(prop as fn(Samples, u8) -> bool);
    }

    #[test]
    fn test_statistic_is_n_div_length_for_sample_with_additional_n_low_and_high_values() {
        fn prop(xs: Samples, n: u8) -> bool {
            // Add extra sample of early weight to ys.
            let min = xs.min();
            let max = xs.max();
            let mut ys = xs.clone();
            for j in 0..n {
                ys.vec.push(min - (j as u64) - 1);
                ys.vec.push(max + (j as u64) + 1);
            }

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = n as f64 / ys.vec.len() as f64;

            (result.statistic - expected).abs() < EPSILON
        }

        check(prop as fn(Samples, u8) -> bool);
    }

    #[test]
    fn test_statistic_is_n_or_m_div_length_for_sample_with_additional_n_low_and_m_high_values() {
        fn prop(xs: Samples, n: u8, m: u8) -> bool {
            // Add extra sample of early weight to ys.
            let min = xs.min();
            let max = xs.max();
            let mut ys = xs.clone();
            for j in 0..n {
                ys.vec.push(min - (j as u64) - 1);
            }

            for j in 0..m {
                ys.vec.push(max + (j as u64) + 1);
            }

            let result = test(&xs.vec, &ys.vec, 0.95);
            let expected = cmp::max(n, m) as f64 / ys.vec.len() as f64;

            (result.statistic - expected).abs() < EPSILON
        }

        check(prop as fn(Samples, u8, u8) -> bool);
    }
}
