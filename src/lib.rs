/// Empirical cumulative distribution function.
pub struct Ecdf<T: Ord> {
  samples: Vec<T>,
  length: usize
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
  pub fn new(samples: &Vec<T>) -> Ecdf<T> {
    let length = samples.len();
    assert!(length > 0);

    // Sort a cloned sample for binary searching.
    let mut sorted = samples.clone();
    sorted.sort();

    Ecdf {
      samples: sorted,
      length: length
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
    let mut index = match self.samples.binary_search(&t) {
      Ok(index) => index,
      Err(index) => index,
    };

    // Walk down samples until at last index that == t.
    while index + 1 < self.length && self.samples[index + 1] == t {
      index += 1;
    }

    // Compensate for 0-based indexing
    if index < self.length && self.samples[index] == t {
      index += 1;
    }

    (index as f64) / (self.length as f64)
  }
}

#[cfg(test)]
extern crate quickcheck;

#[cfg(test)]
mod tests {
  use super::Ecdf;
  use quickcheck::quickcheck;
  use quickcheck::TestResult;

  #[test]
  #[should_panic]
  fn ecdf_panics_on_empty_sample_set() {
    let samples: Vec<i64> = vec!();
    Ecdf::new(&samples);
  }

  #[test]
  fn ecdf_between_zero_and_one() {
    fn prop(samples: Vec<i64>, val: i64) -> TestResult {
      match samples.len() {
        0 => // Discard empty Vec.
          TestResult::discard(),
        _ => TestResult::from_bool({
          let ecdf = Ecdf::new(&samples);
          let actual = ecdf.value(val);

          0.0 <= actual && actual <= 1.0
        })
      }
    }

    quickcheck(prop as fn(Vec<i64>, i64) -> TestResult);
  }

  #[test]
  fn ecdf_is_an_increasing_function() {
    fn prop(samples: Vec<i64>, val: i64) -> TestResult {
      match samples.len() {
        0 => // Discard empty Vec.
          TestResult::discard(),
        _ => TestResult::from_bool({
          let ecdf = Ecdf::new(&samples);
          let actual = ecdf.value(val);

          ecdf.value(val - 1) <= actual && actual <= ecdf.value(val + 1)
        })
      }
    }

    quickcheck(prop as fn(Vec<i64>, i64) -> TestResult);
  }

  #[test]
  fn ecdf_sample_min_minus_one_is_zero() {
    fn prop(samples: Vec<i64>) -> TestResult {
      match samples.iter().min() {
        None => // Discard empty Vec.
          TestResult::discard(),
        Some(&min) => TestResult::from_bool({
          let ecdf = Ecdf::new(&samples);
          ecdf.value(min - 1) == 0.0
        })
      }
    }

    quickcheck(prop as fn(Vec<i64>) -> TestResult);
  }

  #[test]
  fn ecdf_sample_max_is_one() {
    fn prop(samples: Vec<i64>) -> TestResult {
      match samples.iter().max() {
        None => // Discard empty Vec.
          TestResult::discard(),
        Some(&max) => TestResult::from_bool({
          let ecdf = Ecdf::new(&samples);
          ecdf.value(max) == 1.0
        })
      }
    }

    quickcheck(prop as fn(Vec<i64>) -> TestResult);
  }

  #[test]
  fn ecdf_sample_val_is_num_samples_less_than_or_equal_val_div_length() {
    fn prop(samples: Vec<i64>) -> TestResult {
      match samples.first() {
        None => // Discard empty Vec.
          TestResult::discard(),
        Some(&val) => TestResult::from_bool({
          let num_samples_less_than_or_equal_val = samples.iter().filter(|&&x| x <= val).count();
          let expected = (num_samples_less_than_or_equal_val as f64) / (samples.len() as f64);

          let ecdf = Ecdf::new(&samples);

          ecdf.value(val) == expected
        })
      }
    }

    quickcheck(prop as fn(Vec<i64>) -> TestResult);
  }

  #[test]
  fn ecdf_non_sample_val_is_num_samples_less_than_or_equal_val_div_length() {
    fn prop(samples: Vec<i64>, val: i64) -> TestResult {
      let length = samples.len();

      if length == 0 {
        // Discard empty Vec.
        TestResult::discard()
      } else if samples.iter().any(|&x| x == val) {
        // Discard Vec containing val.
        TestResult::discard()
      } else {
        TestResult::from_bool({
          let num_samples_less_than_or_equal_val = samples.iter().filter(|&&x| x <= val).count();
          let expected = (num_samples_less_than_or_equal_val as f64) / (length as f64);

          let ecdf = Ecdf::new(&samples);

          ecdf.value(val) == expected
        })
      }
    }

    quickcheck(prop as fn(Vec<i64>, i64) -> TestResult);
  }
}
