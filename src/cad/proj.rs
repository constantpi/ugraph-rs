use num::BigRational;

use crate::polynomial::{Exponent, Polynomial};

/// Polynomialから、最高次の項の次数と係数を取り出す関数。
fn collect_terms(poly: &Polynomial) -> Option<Vec<(u32, Exponent, BigRational)>> {
    poly.raw_iter()
        .map(|(exp, coeff)| {
            exp.split_last()
                .map(|(last, rest)| (last, rest, coeff.clone()))
        })
        .collect()
}

/// Polynomialから末尾の変数についての係数を取り出す関数
fn collect_last_variable_coefficients(poly: &Polynomial) -> Option<Vec<Polynomial>> {
    let terms = collect_terms(poly)?;
    let max_degree = terms.iter().map(|(last, _, _)| *last).max()?;
    let mut coeffs = vec![Polynomial::zero(); (max_degree + 1) as usize];
    for (last, rest, coeff) in terms {
        let index = last as usize;
        let mut coeff_poly = Polynomial::zero();
        coeff_poly.add_term(rest, coeff);
        coeffs[index] += coeff_poly;
    }

    Some(coeffs)
}

/// 係数の列から微分を計算する関数
fn differentiate_coefficients(coeffs: &[Polynomial]) -> Vec<Polynomial> {
    coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, coeff)| coeff.mul_rational(BigRational::from_integer(i.into())))
        .collect()
}

/// 多項式が末尾の変数について定数であるかどうかを判定する関数。
fn is_constant_for_last_variable(poly: &Polynomial) -> bool {
    if let Some(terms) = collect_terms(poly) {
        terms.iter().all(|(last, _, _)| *last == 0)
    } else {
        true
    }
}
