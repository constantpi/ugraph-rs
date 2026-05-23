use num::{BigRational, One, Zero};
use vec1::{Vec1, vec1};

use super::{UnivariatePolynomial, uni_poly_derivative, uni_poly_div, uni_poly_remainder};
use crate::polyfactor::rational_factorization;

/// 解を表す構造体
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Root {
    poly: UnivariatePolynomial, // 根を求めたい多項式
    range: Range,               // 根の存在範囲
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Range {
    Exact(BigRational),                 // 正確な値がわかっている場合
    Interval(BigRational, BigRational), // 区間で表す場合 (a, b]
}

impl Root {
    pub fn new(
        poly: UnivariatePolynomial,
        lower_bound: BigRational,
        upper_bound: BigRational,
    ) -> Self {
        let lt = poly.leading_coeff();
        // monicである必要がある。
        if !lt.is_one() {
            panic!("Leading coefficient must be 1 or 0 for creating a Root");
        }
        match poly.degree() {
            0 => panic!("The polynomial must have at least one root"),
            1 => {
                let coeffs = poly.iter().cloned().collect::<Vec<_>>();
                let a = coeffs[1].clone();
                let b = coeffs[0].clone();
                if a.is_zero() {
                    panic!("Leading coefficient cannot be zero for a linear polynomial");
                }
                let root = -b / a;
                Root {
                    poly,
                    range: Range::Exact(root),
                }
            }
            _ => {
                if upper_bound <= lower_bound {
                    panic!("Upper bound must be greater than or equal to lower bound");
                }
                Root {
                    poly,
                    range: Range::Interval(lower_bound, upper_bound),
                }
            }
        }
    }

    pub fn new_rational(r: BigRational) -> Self {
        Root {
            poly: UnivariatePolynomial::new(vec1![-r.clone(), BigRational::one()]),
            range: Range::Exact(r),
        }
    }

    pub fn is_same_root(&self, other: &Root) -> bool {
        if self.poly != other.poly {
            false
        } else {
            match (&self.range, &other.range) {
                (Range::Exact(r1), Range::Exact(r2)) => r1 == r2,
                (Range::Interval(l1, u1), Range::Interval(l2, u2)) => {
                    // 区間が重なっているかどうかで同じ根かどうかを判断する
                    l1 < u2 && l2 < u1
                }
                _ => false,
            }
        }
    }

    pub fn get_range(&self) -> &Range {
        &self.range
    }

    pub fn get_poly(&self) -> &UnivariatePolynomial {
        &self.poly
    }

    pub fn get_interval(&self) -> Option<(BigRational, BigRational)> {
        match &self.range {
            Range::Exact(_) => None,
            Range::Interval(lower, upper) => Some((lower.clone(), upper.clone())),
        }
    }

    pub fn get_lower_bound(&self) -> BigRational {
        match &self.range {
            Range::Exact(root) => root.clone(),
            Range::Interval(lower, _) => lower.clone(),
        }
    }

    pub fn get_upper_bound(&self) -> BigRational {
        match &self.range {
            Range::Exact(root) => root.clone(),
            Range::Interval(_, upper) => upper.clone(),
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

fn find_all_roots(poly: &UnivariatePolynomial) -> Vec<Root> {
    let square_free = square_free_part(poly);
    let factors = rational_factorization(&square_free);
    let mut ans = Vec::new();
    for factor in &factors {
        let sequence = eulers_algorithm(factor);
        let (lower_bound, upper_bound, num_of_real_roots) = range_of_roots(&sequence);
        // ここからは根を二分探索で求める
        let mut candidates = vec![(lower_bound, upper_bound, num_of_real_roots)];
        let mut roots = Vec::new();
        let roots = loop {
            let Some((low, high, count)) = candidates.pop() else {
                break roots;
            };
            if count == 1 {
                let root = Root::new(factor.clone(), low, high);
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
        };
        ans.extend(roots);
    }
    ans
}

fn unique_roots(roots: Vec<Root>) -> Vec<Root> {
    let mut unique: Vec<Root> = Vec::new();
    for root in roots {
        if !unique.iter().any(|r| r.is_same_root(&root)) {
            unique.push(root);
        }
    }
    unique
}

/// 一変数の複数の多項式からUniqueなRootをすべて見つける関数
pub fn find_unique_roots(polynomials: &[UnivariatePolynomial]) -> Vec<Root> {
    let all_roots = polynomials
        .iter()
        .filter(|p| p.degree() > 0) // 定数多項式は根を持たないので無視する
        .flat_map(find_all_roots)
        .collect::<Vec<_>>();
    unique_roots(all_roots)
}

impl std::fmt::Display for Root {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.range {
            Range::Exact(root) => write!(f, "{root}"),
            Range::Interval(lower, upper) => write!(f, "[{}, {}] ({})", lower, upper, self.poly),
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
        let roots = find_all_roots(&poly);
        let ans = vec![
            Root::new(
                UnivariatePolynomial::new(vec1![
                    BigRational::from_integer((-1).into()),
                    BigRational::from_integer(1.into())
                ]),
                BigRational::from_integer((-1).into()),
                BigRational::from_integer(1.into()),
            ),
            Root::new(
                UnivariatePolynomial::new(vec1![
                    BigRational::from_integer(1.into()),
                    BigRational::from_integer(1.into())
                ]),
                BigRational::from_integer((-2).into()),
                BigRational::from_integer(2.into()),
            ),
        ];
        assert_eq!(roots, ans);

        // 解が存在しない場合も正しく動くか
        let poly = UnivariatePolynomial::new(vec1![
            BigRational::from_integer(1.into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(1.into())
        ]);
        let roots = find_all_roots(&poly);
        let ans = vec![];
        assert_eq!(roots, ans);
    }
}
