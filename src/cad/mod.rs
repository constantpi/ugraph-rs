mod algorithm;
mod matrix;
mod proj;
mod sturm;
mod univariate;

pub use proj::project_polynomial;
use sturm::{Root, find_all_roots};
pub use univariate::UnivariatePolynomial;
use univariate::{uni_poly_derivative, uni_poly_div, uni_poly_remainder};
