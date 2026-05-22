use itertools::Itertools;
use num::BigInt;
use num_traits::One;
use vec1::Vec1;

use super::BigIntPoly;

/// BigIntPolyの各係数をmod mで-m/2からm/2の範囲に収める関数。中心化処理
fn center_coefficients(f: &BigIntPoly, m: &BigInt) -> BigIntPoly {
    fn center(x: &BigInt, m: &BigInt) -> BigInt {
        let r = x % m;
        let half = m / 2;
        if r > half {
            r - m
        } else if r < -half {
            r + m
        } else {
            r
        }
    }
    let terms = f.get_terms().clone();
    let (rest, last) = terms.split_off_last();
    let rest = rest.iter().map(|c| center(&c, m)).collect::<Vec<_>>();
    let last = center(&last, m);
    let coeffs = Vec1::from_vec_push(rest, last);
    BigIntPoly::new(coeffs)
}

/// 因数分解の候補を組み合わせる関数
pub fn reconstruct_factors(
    f: BigIntPoly,
    factors: Vec1<BigIntPoly>,
    m: &BigInt,
) -> Vec1<BigIntPoly> {
    let n = factors.len();
    for k in 1..=n / 2 {
        // k個の因数の組み合わせを全て試す
        for combination in factors.iter().combinations(k) {
            // 組み合わせに含まれる因数を掛け合わせる
            let product = combination
                .iter()
                .fold(BigIntPoly::one(), |acc, factor| acc * (*factor).clone());
            let centered_product = center_coefficients(&product, m);
            let rem = f.clone() % centered_product.clone();
            if rem.is_zero() {
                let div = f.clone() / centered_product.clone();
                // factorsからcombinationを除いたもの
                let remaining_factors = factors
                    .iter()
                    .filter(|factor| !combination.contains(factor))
                    .cloned()
                    .collect::<Vec<_>>();
                // 再帰的に関数を呼び出す
                let mut result =
                    reconstruct_factors(div, Vec1::try_from_vec(remaining_factors).unwrap(), m);
                result.push(centered_product);
                return result;
            }
        }
    }
    // どの組み合わせもfを割り切らなかった場合は、f自体が素因数であるとみなす
    Vec1::new(f)
}
