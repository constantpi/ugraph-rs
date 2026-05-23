use num::{BigRational, One, Zero};
use vec1::Vec1;

use crate::cad::univariate::{
    UnivariatePolynomial, uni_poly_derivative, uni_poly_div, uni_poly_remainder,
};

/// 解を表す構造体
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Root {
    poly: UnivariatePolynomial, // 根を求めたい多項式
    upper_bound: BigRational,   // 根の上界
    lower_bound: BigRational,   // 根の下界
}

impl Root {
    pub fn new(
        poly: UnivariatePolynomial,
        upper_bound: BigRational,
        lower_bound: BigRational,
    ) -> Self {
        Root {
            poly,
            upper_bound,
            lower_bound,
        }
    }
}

/// オイラーの互除法
fn eulers_algorithm(poly: &UnivariatePolynomial) -> Vec1<UnivariatePolynomial> {
    let mut sequence = Vec1::new(poly.clone());
    let derivative = uni_poly_derivative(poly);
    sequence.push(derivative.clone());
    let (mut f, mut g) = (poly.clone(), derivative);
    loop {
        let r = uni_poly_remainder(&f, &g);
        if let Some(r) = r
            && !r.is_zero()
        {
            sequence.push(-r.clone());
            f = g;
            g = -r;
        } else {
            break sequence;
        }
    }
}

fn square_free_part(poly: &UnivariatePolynomial) -> UnivariatePolynomial {
    let deg = poly.degree();
    if deg <= 1 {
        poly.clone()
    } else {
        let sequence = eulers_algorithm(poly);
        let last = sequence.last();
        uni_poly_div(poly, last).unwrap().monic()
    }
}

/// 列の多項式を受け取って、値を代入し正負の変化を数える関数
fn count_sign_changes(sequence: &[UnivariatePolynomial], value: &BigRational) -> usize {
    let values = sequence
        .iter()
        .map(|p| p.substitute(value))
        .collect::<Vec<_>>();
    let mut sign_changes = 0;
    let mut prev_sign = None;
    for value in values {
        let sign = value.cmp(&BigRational::zero());
        if sign == std::cmp::Ordering::Equal {
            continue; // 0の場合は符号の変化にカウントしない
        }
        if let Some(prev) = prev_sign
            && prev != sign
        {
            sign_changes += 1;
        }
        prev_sign = Some(sign);
    }
    sign_changes
}

/// 列の多項式を受け取って、無限大に近づけたときの符号の変化を数える関数
fn count_sign_changes_at_infinity(
    sequence: &[UnivariatePolynomial],
    is_positive_infinity: bool,
) -> usize {
    let leading_coeffs = sequence
        .iter()
        .map(|p| p.leading_coeff())
        .collect::<Vec<_>>();
    let mut sign_changes = 0;
    let mut prev_sign = None;
    let mut reversed = false;
    for coeff in leading_coeffs {
        let sign = coeff.cmp(&BigRational::zero());
        if sign == std::cmp::Ordering::Equal {
            panic!("Leading coefficient cannot be zero for counting sign changes at infinity");
        }
        let sign = if reversed { sign.reverse() } else { sign };
        if let Some(prev) = prev_sign
            && prev != sign
        {
            sign_changes += 1;
        }
        prev_sign = Some(sign);
        if !is_positive_infinity {
            reversed = !reversed;
        }
    }
    sign_changes
}

/// 実数解の個数を求める関数
fn count_real_roots(sequence: &[UnivariatePolynomial]) -> usize {
    let sign_changes_at_positive_infinity = count_sign_changes_at_infinity(sequence, true);
    let sign_changes_at_negative_infinity = count_sign_changes_at_infinity(sequence, false);
    sign_changes_at_negative_infinity - sign_changes_at_positive_infinity
}

/// 指定された範囲の根の個数を数える関数
fn count_roots_in_range(
    sequence: &[UnivariatePolynomial],
    low: &BigRational,
    high: &BigRational,
) -> usize {
    let sign_changes_low = count_sign_changes(sequence, low);
    let sign_changes_high = count_sign_changes(sequence, high);
    sign_changes_low - sign_changes_high
}

/// 解の存在範囲を求める関数
fn range_of_roots(sequence: &Vec1<UnivariatePolynomial>) -> (BigRational, BigRational, usize) {
    let number_of_real_roots = count_real_roots(sequence);
    let mut abs = BigRational::one();
    loop {
        let low = -&abs;
        if number_of_real_roots == count_roots_in_range(sequence, &low, &abs) {
            break (low, abs, number_of_real_roots);
        }
        abs *= BigRational::from_integer(2.into());
    }
}

pub fn sturm(poly: &UnivariatePolynomial) -> Vec<Root> {
    let square_free = square_free_part(poly);
    let sequence = eulers_algorithm(&square_free);
    let (lower_bound, upper_bound, num_of_real_roots) = range_of_roots(&sequence);
    // ここからは根を二分探索で求める
    let mut candidates = vec![(lower_bound, upper_bound, num_of_real_roots)];
    let mut roots = Vec::new();
    loop {
        let Some((low, high, count)) = candidates.pop() else {
            break roots;
        };
        if count == 1 {
            let root = Root::new(poly.clone(), high, low);
            roots.push(root);
            continue;
        }

        let mid = (&low + &high) / BigRational::from_integer(2.into());
        let low_count = count_roots_in_range(&sequence, &low, &mid);
        let high_count = count - low_count;
        if low_count > 0 {
            candidates.push((low, mid.clone(), low_count));
        }
        if high_count > 0 {
            candidates.push((mid, high, high_count));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_square_free_part() {
        // (x+1)(x+1)(x-1) = x^3 + x^2 - x - 1 の平方因子を取り除くと x^2 - 1 になるはず
        let poly = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-1).into()),
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(1.into()),
            BigRational::from_integer(1.into())
        ]);
        let result = square_free_part(&poly);
        let expected = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(1.into())
        ]);
        assert_eq!(result, expected);
    }

    #[test]
    fn test_sturm() {
        // x^3 + x^2 - x - 1 の根は -1 と 1 なので、Sturmの定理で正しく求まるか
        let poly = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-1).into()),
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(1.into()),
            BigRational::from_integer(1.into())
        ]);
        let roots = sturm(&poly);
        let ans = vec![
            Root::new(
                poly.clone(),
                BigRational::from_integer(2.into()),
                BigRational::from_integer(0.into()),
            ),
            Root::new(
                poly.clone(),
                BigRational::from_integer(0.into()),
                BigRational::from_integer((-2).into()),
            ),
        ];
        assert_eq!(roots, ans);

        // 解が存在しない場合も正しく動くか
        let poly = UnivariatePolynomial::new(vec1![
            BigRational::from_integer(1.into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(1.into())
        ]);
        let roots = sturm(&poly);
        let ans = vec![];
        assert_eq!(roots, ans);
    }
}
