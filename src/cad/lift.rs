use color_eyre::{Result, eyre::eyre};

use super::{Root, find_unique_roots, specialize_polynomial};
use crate::polynomial::Polynomial;

/// あるサンプル点に対して多変数多項式の解を見つける関数
fn find_solutions_at_sample_point(poly_list: &[Polynomial], sample: &[Root]) -> Result<Vec<Root>> {
    let values = sample
        .iter()
        .map(|root| root.get_poly())
        .cloned()
        .collect::<Vec<_>>();
    let specialized_polys = poly_list
        .iter()
        .map(|poly| specialize_polynomial(poly, &values))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| eyre!("Failed to specialize polynomials at the given sample points"))?;
    Ok(find_unique_roots(&specialized_polys))
}

/// CADのLiftステップを実装する関数
/// あるサンプル点に対して、そこから上の次元の解を見つける関数
pub fn lifting(poly_list: &[Polynomial], sample_list: &[Vec<Root>]) -> Result<Vec<Vec<Root>>> {
    sample_list
        .iter()
        .map(|sample| find_solutions_at_sample_point(poly_list, sample))
        .collect::<Result<Vec<_>>>()
}
