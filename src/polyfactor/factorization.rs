use num::{BigInt, BigRational, Integer};
use num_traits::{One, Zero};
use vec1::Vec1;

use super::{BigIntPoly, berlekamp_factorization, hensel_lifting, reconstruct_factors};
use crate::cad::UnivariatePolynomial;

pub fn rational_factorization(poly: &UnivariatePolynomial) -> Vec1<UnivariatePolynomial> {
    let coeffs = poly.get_coeffs();
    let integer_coeffs = rational_to_integer_coeffs(coeffs);
    let (monic_coeffs, leading_coeff) = monicize(&integer_coeffs);
    let poly = BigIntPoly::new(monic_coeffs.clone());
    let prime_factors = berlekamp_factorization(monic_coeffs);
    let (lifted_factors, modulo) = hensel_lifting(&poly, &prime_factors);
    let factors = reconstruct_factors(poly, lifted_factors, &modulo);
    let (rest, last) = factors.split_off_last();
    let rest = rest
        .iter()
        .map(|factor| demonicize_to_univariate(factor, &leading_coeff))
        .collect::<Vec<_>>();
    let last = demonicize_to_univariate(&last, &leading_coeff);
    Vec1::from_vec_push(rest, last)
}

/// BigRationalの係数からBigIntの係数に変換する関数
fn rational_to_integer_coeffs(coeffs: &[BigRational]) -> Vec<BigInt> {
    // まずは分母の最小公倍数を求める
    let lcm = coeffs.iter().fold(BigInt::from(1), |acc, coeff| {
        let denom = coeff.denom();
        acc.lcm(denom)
    });
    // 各係数を最小公倍数で割って整数に変換
    coeffs
        .iter()
        .map(|coeff| {
            let numerator = coeff.numer();
            numerator * &lcm / coeff.denom()
        })
        .collect()
}

/// 多項式のモニック化
fn monicize(coeffs: &[BigInt]) -> (Vec1<BigInt>, BigInt) {
    if let Some(leading) = coeffs.last() {
        // a^(n-1) * f(x/a) を計算する
        if leading.is_zero() {
            panic!("Leading coefficient cannot be zero for monicization");
        }
        let n = coeffs.len() - 1;
        let rest = coeffs
            .iter()
            .enumerate()
            .take(n)
            .map(|(i, coeff)| coeff * leading.pow((n - 1 - i) as u32))
            .collect();
        (Vec1::from_vec_push(rest, BigInt::one()), leading.clone())
    } else {
        (Vec1::new(BigInt::from(0)), BigInt::one())
    }
}

/// モニック化したときに変数変換をしているのでそれを元に戻す関数
fn demonicize(poly: &BigIntPoly, leading_coeff: &BigInt) -> Vec1<BigInt> {
    let terms = poly.get_terms().clone();
    let (first, rest) = terms.split_off_first();
    let rest = rest
        .iter()
        .enumerate()
        .map(|(i, coeff)| coeff * leading_coeff.pow((i + 1) as u32))
        .collect::<Vec<_>>();
    let mut coeff = Vec1::new(first);
    coeff.extend(rest);
    coeff
}

/// BigIntPolyをUnivariatePolynomialに変換する関数
/// その際にモニックにする
fn bigintpoly_to_univariate(poly: &Vec1<BigInt>) -> UnivariatePolynomial {
    let (rest, lt) = poly.clone().split_off_last();
    let rest = rest
        .iter()
        .map(|coeff| BigRational::new(coeff.clone(), lt.clone()))
        .collect::<Vec<_>>();
    let coeffs = Vec1::from_vec_push(rest, BigRational::one());
    UnivariatePolynomial::new(coeffs)
}

/// demonicizeした後の多項式をUnivariatePolynomialに変換する関数
fn demonicize_to_univariate(poly: &BigIntPoly, leading_coeff: &BigInt) -> UnivariatePolynomial {
    let demonicized_coeffs = demonicize(poly, leading_coeff);
    bigintpoly_to_univariate(&demonicized_coeffs)
}

#[cfg(test)]
mod tests {
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_rational_to_integer_coeffs() {
        let coeffs = vec![
            BigRational::new(BigInt::from(1), BigInt::from(2)),
            BigRational::new(BigInt::from(5), BigInt::from(6)),
            BigRational::new(BigInt::from(7), BigInt::from(4)),
        ];
        let expected = vec![BigInt::from(6), BigInt::from(10), BigInt::from(21)];
        assert_eq!(rational_to_integer_coeffs(&coeffs), expected);
    }

    #[test]
    fn test_monicize() {
        let coeffs = vec![BigInt::from(2), BigInt::from(3), BigInt::from(4)];
        let expected = vec1![BigInt::from(8), BigInt::from(3), BigInt::from(1)];
        assert_eq!(monicize(&coeffs), (expected, BigInt::from(4)));
    }

    #[test]
    fn test_rational_factorization() {
        let coeffs = vec1![
            BigRational::new(BigInt::from(-3), BigInt::from(5)),
            BigRational::new(BigInt::from(-6), BigInt::from(5)),
            BigRational::new(BigInt::from(1), BigInt::from(10)),
            BigRational::new(BigInt::from(1), BigInt::from(5)),
        ];
        let poly = UnivariatePolynomial::new(coeffs);
        let factors = rational_factorization(&poly);
        let expected = vec1![
            UnivariatePolynomial::new(vec1![
                BigRational::new(BigInt::from(-6), BigInt::from(1)),
                BigRational::new(BigInt::from(0), BigInt::from(1)),
                BigRational::new(BigInt::from(1), BigInt::from(1))
            ]),
            UnivariatePolynomial::new(vec1![
                BigRational::new(BigInt::from(1), BigInt::from(2)),
                BigRational::new(BigInt::from(1), BigInt::from(1))
            ])
        ];
        assert_eq!(factors, expected);
    }
}
