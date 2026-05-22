use std::vec;

use num::{BigInt, BigRational, Integer};
use num_traits::{One, Zero};
use vec1::Vec1;

use super::{BigIntPoly, berlekamp_factorization, hensel_lifting};

pub fn rational_factorization(coeffs: Vec<BigRational>) -> Vec<BigIntPoly> {
    let integer_coeffs = rational_to_integer_coeffs(&coeffs);
    let monic_coeffs = monicize(&integer_coeffs);
    let poly = BigIntPoly::new(monic_coeffs.clone());
    // println!("Monic polynomial: {:?}", poly);
    let prime_factors = berlekamp_factorization(monic_coeffs);
    let (lifted_factors, modulo) = hensel_lifting(&poly, &prime_factors);
    // println!("modulo: {}", modulo);
    // for factor in lifted_factors {
    //     println!("Factor: {:?}", factor);
    // }

    vec![]
}

/// BigRationalの係数からBigIntの係数に変換する関数
fn rational_to_integer_coeffs(coeffs: &Vec<BigRational>) -> Vec<BigInt> {
    // まずは分母の最小公倍数を求める
    let lcm = coeffs.iter().fold(BigInt::from(1), |acc, coeff| {
        let denom = coeff.denom();
        acc.lcm(denom)
    });
    // 各係数を最小公倍数で割って整数に変換
    coeffs
        .into_iter()
        .map(|coeff| {
            let numerator = coeff.numer();
            numerator * &lcm / coeff.denom()
        })
        .collect()
}

/// 多項式のモニック化
fn monicize(coeffs: &Vec<BigInt>) -> Vec1<BigInt> {
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
        Vec1::from_vec_push(rest, BigInt::one())
    } else {
        Vec1::new(BigInt::from(0))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
        let expected = vec![BigInt::from(8), BigInt::from(3), BigInt::from(1)];
        assert_eq!(monicize(&coeffs), expected);
    }

    #[test]
    fn test_rational_factorization() {
        let coeffs = vec![
            BigRational::new(BigInt::from(-3), BigInt::from(5)),
            BigRational::new(BigInt::from(-6), BigInt::from(5)),
            BigRational::new(BigInt::from(1), BigInt::from(10)),
            BigRational::new(BigInt::from(1), BigInt::from(5)),
        ];
        let factors = rational_factorization(coeffs);
    }
}
