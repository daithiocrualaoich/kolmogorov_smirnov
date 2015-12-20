extern crate kolmogorov_smirnov as ks;

use ks::calculate_critical_value;

use std::env;

/// Calculate critical values dataset.
///
/// # Examples
///
/// ```bash
/// cargo run --bin critical_values <confidence> <num_samples> <limit>
/// ```
///
/// This will print the critical values of the Kolmogorov-Smirnov two sample
/// test for samples of size `<num_samples>` against samples of sizes 16
/// through `<limit>` inclusive at the specified confidence level.
///
/// `<num_samples>` and `<limit>` must be positive integers, `<confidence>` must
/// be a floating point number strictly between zero and one.
fn main() {
    let args: Vec<String> = env::args().collect();

    let confidence: f64 = args[1].parse().expect("<confidence> must be a floating point number.");
    let n1: usize = args[2].parse().expect("<num_samples> must be an integer.");
    let limit: usize = args[3].parse().expect("<limit> must be an integer.");

    assert!(n1 > 0 && limit > 0);
    assert!(0.0 < confidence && confidence < 1.0);

    println!("n1\tn2\tconfidence\tcritical_value");
    for n2 in 16..(limit + 1) {
        println!("{}\t{}\t{}\t{}",
                 n1,
                 n2,
                 confidence,
                 calculate_critical_value(n1, n2, confidence));
    }
}
