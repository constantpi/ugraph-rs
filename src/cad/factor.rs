use num::{BigRational, Zero};
use vec1::Vec1;

use super::UnivariatePolynomial;
use crate::polyfactor::rational_factorization;
use crate::polynomial::{Exponent, Polynomial};

/// もし一変数の多項式であれば因数分解を行う関数
pub fn factorization(poly: &Polynomial) -> Option<Vec<Polynomial>> {
    // まずは一変数の多項式かどうかを確認する
    let variables = poly
        .raw_iter()
        .map(|(exps, _)| {
            exps.iter()
                .enumerate()
                .filter(|&(_, exp)| exp > 0)
                .map(|(i, _)| i)
        })
        .flatten()
        .collect::<Vec<_>>();
    let Some(var) = variables.first() else {
        // 変数がない場合は定数多項式なので、因数分解はできない
        return None;
    };
    if variables.iter().any(|&other| other != *var) {
        // 変数が複数ある場合は因数分解はできない
        return None;
    }
    let total_num_vars = poly.total_num_vars();
    let mut coeffs = vec![];
    for (exps, coeff) in poly.raw_iter() {
        let exps = exps.iter().collect::<Vec<_>>();
        if let Some(exp) = exps.get(*var as usize) {
            // coeffs[exp] = coeffとしたいが、expは0から始まるとは限らないので、必要に応じてcoeffsを拡張する
            while coeffs.len() <= *exp as usize {
                coeffs.push(BigRational::zero());
            }
            coeffs[*exp as usize] = coeff.clone();
        }
    }
    let Ok(coeffs) = Vec1::try_from_vec(coeffs) else {
        return None;
    };
    let uni_poly = UnivariatePolynomial::new(coeffs);
    let factors = rational_factorization(&uni_poly);
    let factors_poly = factors
        .into_iter()
        .map(|factor| {
            let mut polyfactor = Polynomial::zero();
            for (i, coeff) in factor.iter().enumerate() {
                let mut exps = vec![0; total_num_vars];
                exps[*var as usize] = i as u32;
                polyfactor.add_term(Exponent::new(exps), coeff.clone());
            }
            polyfactor
        })
        .collect();
    Some(factors_poly)
}
