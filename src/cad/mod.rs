mod algorithm;
mod binary_search;
mod lift;
mod matrix;
mod proj;
mod specialize;
mod sturm;
mod univariate;

pub use algorithm::find_solution;
pub use binary_search::{calc_sample_points, refine_root};
use lift::lifting;
pub use proj::project_polynomial;
use proj::psc_0;
use specialize::specialize_polynomial;
use sturm::{Range, Root, find_unique_roots};
pub use univariate::{UnivariatePolynomial, polynomial_to_univariate};
use univariate::{uni_poly_derivative, uni_poly_div, uni_poly_remainder};
