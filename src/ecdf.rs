//! Empirical cumulative distribution function.

pub struct Ecdf<T: Ord> {
    samples: Vec<T>,
    length: usize,
}

impl<T: Ord + Clone> Ecdf<T> {
    /// Construct a new representation of a cumulative distribution function for
    /// a given sample.
    ///
    /// The construction will involve computing a sorted clone of the given sample
    /// and may be inefficient or completely prohibitive for large samples. This
    /// computation is amortized significantly if there is heavy use of the value
    /// function.
    ///
    /// # Panics
    ///
    /// The sample set must be non-empty.
    ///
    /// # Examples
    ///
    /// ```
    /// extern crate kolmogorov_smirnov as ks;
    ///
    /// let samples = vec!(9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    /// let ecdf = ks::Ecdf::new(&samples);
    /// ```
    pub fn new(samples: &[T]) -> Ecdf<T> {
        let length = samples.len();
        assert!(length > 0);

        // Sort a copied sample for binary searching.
        let mut sorted = samples.to_vec();
        sorted.sort();

        Ecdf {
            samples: sorted,
            length: length,
        }
    }

    /// Calculate a value of the empirical cumulative distribution function for
    /// a given sample.
    ///
    /// # Examples
    ///
    /// ```
    /// extern crate kolmogorov_smirnov as ks;
    ///
    /// let samples = vec!(9, 8, 7, 6, 5, 4, 3, 2, 1, 0);
    /// let ecdf = ks::Ecdf::new(&samples);
    /// println!("{}", ecdf.value(4));
    /// ```
    pub fn value(&self, t: T) -> f64 {
        let num_samples_leq_t = match self.samples.binary_search(&t) {
            Ok(mut index) => {
                // At least one sample is a t and we have the index of it. Need
                // to walk down the sorted samples until at last that == t.
                while index + 1 < self.length && self.samples[index + 1] == t {
                    index += 1;
                }

                // Compensate for 0-based indexing.
                index + 1
            }
            Err(index) => {
                // No sample is a t but if we had to put one in it would go at
                // index. This means all indices to the left have samples < t
                // and should be counted in the cdf proportion. We must take one
                // from index to get the last included sample but then we just
                // have to add one again to account for 0-based indexing.
                index
            }
        };

        num_samples_leq_t as f64 / self.length as f64
    }
}

#[cfg(test)]
mod tests {
    extern crate quickcheck;
    extern crate rand;

    use self::quickcheck::{quickcheck, TestResult, Arbitrary, Gen};
    use self::rand::Rng;
    use super::Ecdf;

    /// Wrapper for generating sample data with QuickCheck.
    ///
    /// Samples must be non-empty sequences of i64 values.
    #[derive(Debug, Clone)]
    struct Samples {
        vec: Vec<i64>,
    }

    impl Arbitrary for Samples {
        fn arbitrary<G: Gen>(g: &mut G) -> Samples {
            let max = g.size();
            let size = g.gen_range(1, max);
            let vec = (0..size).map(|_| Arbitrary::arbitrary(g)).collect();

            Samples { vec: vec }
        }
    }

    #[test]
    #[should_panic(expected="assertion failed: length > 0")]
    fn ecdf_panics_on_empty_samples_set() {
        let xs: Vec<i64> = vec![];
        Ecdf::new(&xs);
    }

    #[test]
    fn ecdf_between_zero_and_one() {
        fn prop(xs: Samples, val: i64) -> bool {
            let ecdf = Ecdf::new(&xs.vec);
            let actual = ecdf.value(val);

            0.0 <= actual && actual <= 1.0
        }

        quickcheck(prop as fn(Samples, i64) -> bool);
    }

    #[test]
    fn ecdf_is_an_increasing_function() {
        fn prop(xs: Samples, val: i64) -> bool {
            let ecdf = Ecdf::new(&xs.vec);
            let actual = ecdf.value(val);

            ecdf.value(val - 1) <= actual && actual <= ecdf.value(val + 1)
        }

        quickcheck(prop as fn(Samples, i64) -> bool);
    }

    #[test]
    fn ecdf_sample_min_minus_one_is_zero() {
        fn prop(xs: Samples) -> bool {
            let &min = xs.vec.iter().min().unwrap();
            let ecdf = Ecdf::new(&xs.vec);

            ecdf.value(min - 1) == 0.0
        }

        quickcheck(prop as fn(Samples) -> bool);
    }

    #[test]
    fn ecdf_sample_max_is_one() {
        fn prop(xs: Samples) -> bool {
            let &max = xs.vec.iter().max().unwrap();
            let ecdf = Ecdf::new(&xs.vec);

            ecdf.value(max) == 1.0
        }

        quickcheck(prop as fn(Samples) -> bool);
    }

    #[test]
    fn ecdf_sample_val_is_num_samples_leq_val_div_length() {
        fn prop(xs: Samples) -> bool {
            let &val = xs.vec.first().unwrap();
            let num_samples = xs.vec
                                .iter()
                                .filter(|&&x| x <= val)
                                .count();
            let expected = num_samples as f64 / xs.vec.len() as f64;

            let ecdf = Ecdf::new(&xs.vec);

            ecdf.value(val) == expected
        }

        quickcheck(prop as fn(Samples) -> bool);
    }

    #[test]
    fn ecdf_non_sample_val_is_num_samples_leq_val_div_length() {
        fn prop(xs: Samples, val: i64) -> TestResult {
            let length = xs.vec.len();

            if xs.vec.iter().any(|&x| x == val) {
                // Discard Vec containing val.
                return TestResult::discard();
            }

            let num_samples = xs.vec
                                .iter()
                                .filter(|&&x| x <= val)
                                .count();
            let expected = num_samples as f64 / length as f64;

            let ecdf = Ecdf::new(&xs.vec);

            TestResult::from_bool(ecdf.value(val) == expected)
        }

        quickcheck(prop as fn(Samples, i64) -> TestResult);
    }
}
