use num::BigRational;

use crate::cad::{Root, UnivariatePolynomial, polynomial_to_univariate, psc_0};
use crate::polynomial::{Exponent, Polynomial};

fn generate_constant_polynomial(num_vars: usize, constant: &BigRational) -> Polynomial {
    let mut p = Polynomial::zero();
    p.add_term(Exponent::new(vec![0; num_vars]), constant.clone());
    p
}

/// 他変数多項式から1変数目についての多項式を作る関数
fn collect_terms(poly: &Polynomial) -> Option<Vec<(u32, Exponent, BigRational)>> {
    poly.raw_iter()
        .map(|(exp, coeff)| {
            exp.split_first()
                .map(|(first, rest)| (first, rest, coeff.clone()))
        })
        .collect()
}

/// Polynomialから最初の変数についての係数を取り出す関数
/// 例えば、x^2*y + 3*x*y^2 + 5を渡すと、[5, x^2, 3x]を返す
fn collect_first_variable_coefficients(poly: &Polynomial) -> Option<Vec<Polynomial>> {
    let terms = collect_terms(poly)?;
    let max_degree = terms.iter().map(|(last, _, _)| *last).max()?;
    let mut coeffs = vec![Polynomial::zero(); (max_degree + 1) as usize];
    for (last, rest, coeff) in terms {
        let index = last as usize;
        coeffs[index].add_term(rest, coeff);
    }

    Some(coeffs)
}

/// UnivariatePolynomialから係数列を取り出す関数
fn univariate_to_coefficients(poly: &UnivariatePolynomial, num_vars: usize) -> Vec<Polynomial> {
    poly.iter()
        .map(|coeff| generate_constant_polynomial(num_vars, coeff))
        .collect()
}

/// 他変数関数の一番初めの変数に値を代入する関数
fn substitute_first_variable(
    poly: &Polynomial,
    value: &UnivariatePolynomial,
    num_vars: usize, // 残った変数の数
) -> Option<Polynomial> {
    let poly_coeffs = collect_first_variable_coefficients(poly)?;
    let value_coeffs = univariate_to_coefficients(value, num_vars);

    let resultant = psc_0(&poly_coeffs, &value_coeffs);
    Some(resultant)
}

/// 他変数関数に順番に値を代入していく関数
/// 3変数関数ならvaluesは2つのUnivariatePolynomialを持つべき
pub fn specialize_polynomial(
    poly: &Polynomial,
    values: &[UnivariatePolynomial],
) -> Option<UnivariatePolynomial> {
    let current_num_vars = values.len();
    let substituted_poly =
        values
            .iter()
            .enumerate()
            .fold(Some(poly.clone()), |current_poly, (i, value)| {
                current_poly
                    .and_then(|p| substitute_first_variable(&p, value, current_num_vars - i))
            })?;
    polynomial_to_univariate(&substituted_poly).ok()
}

#[cfg(test)]
mod tests {
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_specialize_polynomial() {
        let poly = {
            let mut p = Polynomial::zero();
            p.add_term(
                Exponent::new(vec![2, 2, 2]),
                BigRational::from_integer(1.into()),
            ); // x^2*y^2
            p.add_term(
                Exponent::new(vec![1, 1, 1]),
                BigRational::from_integer(10.into()),
            ); // 10*x*y
            p.add_term(
                Exponent::new(vec![0, 0, 0]),
                BigRational::from_integer(5.into()),
            ); // 5
            p
        };
        let x = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-2).into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(1.into())
        ]); // x = sqrt(2)
        let y = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(1.into()),
        ]); // y = 1
        let substituted = specialize_polynomial(&poly, &[x, y]).unwrap();
        let ans = UnivariatePolynomial::new(vec1![
            BigRational::from_integer(25.into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer((-180).into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(4.into())
        ]); // 4x^4 - 180x^2 + 25
        assert_eq!(substituted, ans);
    }
}
