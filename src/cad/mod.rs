mod algorithm;
mod matrix;
mod proj;
mod sturm;
mod univariate;

pub use algorithm::find_solution;
pub use proj::project_polynomial;
use sturm::{Root, find_all_roots};
pub use univariate::{UnivariatePolynomial, polynomial_to_univariate};
use univariate::{uni_poly_derivative, uni_poly_div, uni_poly_remainder};
