pub mod ecdf;
pub mod test;

pub use test::{test, test_f64, calculate_critical_value};
pub use ecdf::{Ecdf, ecdf, percentile};
