use num::{BigInt, BigRational, Signed, Zero};

use super::{Range, Root, UnivariatePolynomial, refine_range};
use crate::parser::RelOp;
use crate::polyfactor::rational_to_integer_coeffs;
use crate::polynomial::Polynomial;

/// 符号の確定した区間を表す構造体
pub struct SignedInterval {
    is_positive: bool,
    abs_lower: BigRational,
    abs_upper: BigRational,
}

/// 符号の確定した区間または厳密値を表す構造体
pub enum SignedRange {
    Exact(BigRational),
    Interval(SignedInterval),
}

impl SignedInterval {
    /// 区間が正の数のみを含むか
    pub fn is_positive(&self) -> bool {
        self.is_positive
    }

    pub fn get_lower(&self) -> &BigRational {
        &self.abs_lower
    }

    pub fn get_upper(&self) -> &BigRational {
        &self.abs_upper
    }
}

/// lowerとupperが符号をまたがないかどうかを判定する関数
fn is_sign_consistent(lower: &BigRational, upper: &BigRational) -> Option<bool> {
    if lower.is_positive() && upper.is_positive() {
        Some(true)
    } else if lower.is_negative() && upper.is_negative() {
        Some(false)
    } else {
        None
    }
}

/// Rootから符号の確定した区間を取り出す関数
pub fn get_signed_range(root: &Root) -> SignedRange {
    match root.get_range() {
        Range::Exact(r) => SignedRange::Exact(r.clone()),
        Range::Interval(lower, upper) => {
            let poly = root.get_poly();
            let mut lower = lower.clone();
            let mut upper = upper.clone();
            let (lower, upper, is_positive) = loop {
                if let Some(is_positive) = is_sign_consistent(&lower, &upper) {
                    break (lower, upper, is_positive);
                }
                let (new_low, new_high) = refine_range(poly, &lower, &upper, false);
                lower = new_low;
                upper = new_high;
            };
            let abs_lower = lower.abs().min(upper.abs());
            let abs_upper = lower.abs().max(upper.abs());
            SignedRange::Interval(SignedInterval {
                is_positive,
                abs_lower,
                abs_upper,
            })
        }
    }
}

pub fn pow_range(base: &SignedRange, exp: u32) -> SignedRange {
    match base {
        SignedRange::Exact(r) => SignedRange::Exact(r.pow(exp as i32)),
        SignedRange::Interval(interval) => {
            let abs_lower = interval.abs_lower.pow(exp as i32);
            let abs_upper = interval.abs_upper.pow(exp as i32);
            let is_positive = if exp.is_multiple_of(2) {
                true
            } else {
                interval.is_positive
            };
            SignedRange::Interval(SignedInterval {
                is_positive,
                abs_lower,
                abs_upper,
            })
        }
    }
}

fn mul_constant_range(range: &SignedRange, constant: &BigRational) -> SignedRange {
    match range {
        SignedRange::Exact(r) => SignedRange::Exact(r * constant),
        SignedRange::Interval(interval) => {
            let is_positive = if constant.is_positive() {
                interval.is_positive
            } else if constant.is_negative() {
                !interval.is_positive
            } else {
                return SignedRange::Exact(BigRational::zero());
            };
            let abs_lower = interval.abs_lower.clone() * constant.abs();
            let abs_upper = interval.abs_upper.clone() * constant.abs();
            SignedRange::Interval(SignedInterval {
                is_positive,
                abs_lower,
                abs_upper,
            })
        }
    }
}

pub fn mul_ranges(range1: &SignedRange, range2: &SignedRange) -> SignedRange {
    match (range1, range2) {
        (SignedRange::Exact(r1), SignedRange::Exact(r2)) => SignedRange::Exact(r1 * r2),
        (SignedRange::Exact(r), interval) | (interval, SignedRange::Exact(r)) => {
            mul_constant_range(interval, r)
        }
        (SignedRange::Interval(interval1), SignedRange::Interval(interval2)) => {
            let abs_lower = interval1.abs_lower.clone() * interval2.abs_lower.clone();
            let abs_upper = interval1.abs_upper.clone() * interval2.abs_upper.clone();
            let is_positive = interval1.is_positive == interval2.is_positive;
            SignedRange::Interval(SignedInterval {
                is_positive,
                abs_lower,
                abs_upper,
            })
        }
    }
}

/// 最小多項式に対してLandauの不等式によってMahler measureの上界を求める関数
pub fn log_landau_bound(poly: &UnivariatePolynomial) -> usize {
    let coeffs = poly.get_coeffs();
    let integer_coeffs = rational_to_integer_coeffs(coeffs);
    // 係数の二乗和の平方根のlog2を上から取る
    let sum_of_squares = integer_coeffs.iter().map(|c| c * c).sum::<BigInt>();
    let sqrt = sum_of_squares.sqrt() + BigInt::from(1);
    log2_ceil(&sqrt)
}

/// BigIntのlog2を上から取る関数
fn log2_ceil(n: &BigInt) -> usize {
    let mut bound = 0;
    let mut power_of_two = BigInt::from(1);
    loop {
        if power_of_two > *n {
            break bound;
        }
        bound += 1;
        power_of_two *= 2;
    }
}

/// 定数のlog_landau_boundを求める関数
pub fn log_landau_bound_constant(constant: &BigRational) -> Option<usize> {
    if constant.is_zero() {
        None
    } else {
        // 分母と分子のmax
        let max_coeff = std::cmp::max(constant.numer().abs(), constant.denom().abs());
        Some(log2_ceil(&max_coeff))
    }
}

/// Mahler measureと区間
struct SignedMahler {
    dimension: usize,
    log_landau_bound: usize,
    range: SignedRange,
}

/// Rootから符号の確定した区間とLandauの不等式によるMahler measureの上界を求める関数
/// ただし0の場合はNoneを返す
fn get_signed_mahler(root: &Root) -> Option<SignedMahler> {
    let range = get_signed_range(root);
    let dimension = root.get_poly().degree();
    let log_landau_bound = match &range {
        SignedRange::Exact(r) => log_landau_bound_constant(r)?,
        SignedRange::Interval(_) => log_landau_bound(root.get_poly()),
    };

    Some(SignedMahler {
        dimension,
        log_landau_bound,
        range,
    })
}

fn rational_to_signed_mahler(r: &BigRational) -> Option<SignedMahler> {
    let log_landau_bound = log_landau_bound_constant(r)?;
    Some(SignedMahler {
        dimension: 1,
        log_landau_bound,
        range: SignedRange::Exact(r.clone()),
    })
}

fn pow_signed_mahler(mahler: &SignedMahler, exp: u32) -> SignedMahler {
    let range = pow_range(&mahler.range, exp);
    // 累乗しても次元は変わらない
    // boundはexp乗になるので、log_landau_boundはexp倍になる
    SignedMahler {
        dimension: mahler.dimension,
        log_landau_bound: mahler.log_landau_bound * (exp as usize),
        range,
    }
}

fn mul_signed_mahler(mahler1: &SignedMahler, mahler2: &SignedMahler) -> SignedMahler {
    let range = mul_ranges(&mahler1.range, &mahler2.range);
    // 次元は足し算になる
    let dimension = mahler1.dimension * mahler2.dimension;
    // boundはA^d2 * B^d1になるので、log_landau_boundはd2 * log_landau_bound1 + d1 * log_landau_bound2になる
    let log_landau_bound =
        mahler1.log_landau_bound * mahler2.dimension + mahler2.log_landau_bound * mahler1.dimension;
    SignedMahler {
        dimension,
        log_landau_bound,
        range,
    }
}

struct UnsignedMahler {
    lower: BigRational,
    upper: BigRational,
    dimension: usize,
    log_landau_bound: usize,
}

/// SignedMahlerをUnsignedMahlerに変換する関数
fn signed_to_unsigned_mahler(signed: &SignedMahler) -> UnsignedMahler {
    match &signed.range {
        SignedRange::Exact(r) => UnsignedMahler {
            lower: r.clone(),
            upper: r.clone(),
            dimension: signed.dimension,
            log_landau_bound: signed.log_landau_bound,
        },
        SignedRange::Interval(interval) => {
            let is_positive = interval.is_positive();
            let lower = if is_positive {
                interval.get_lower().clone()
            } else {
                -interval.get_upper().clone()
            };
            let upper = if is_positive {
                interval.get_upper().clone()
            } else {
                -interval.get_lower().clone()
            };
            UnsignedMahler {
                lower,
                upper,
                dimension: signed.dimension,
                log_landau_bound: signed.log_landau_bound,
            }
        }
    }
}

fn add_unsigned_mahler(m1: &UnsignedMahler, m2: &UnsignedMahler) -> UnsignedMahler {
    let dimension = m1.dimension + m2.dimension;
    let log_landau_bound =
        m1.log_landau_bound * m2.dimension + m2.log_landau_bound * m1.dimension + dimension;
    UnsignedMahler {
        lower: m1.lower.clone() + m2.lower.clone(),
        upper: m1.upper.clone() + m2.upper.clone(),
        dimension,
        log_landau_bound,
    }
}

pub enum MahlerResult {
    Positive,
    Negative,
    Zero,
    Uncertain,
}

pub fn evaluate_polynomial_by_mahler(poly: &Polynomial, sample: &[Root]) -> MahlerResult {
    let mut term_mahlers = Vec::new();
    for (exp, coeff) in poly.lex_iter() {
        let Some(coeff_mahler) = rational_to_signed_mahler(coeff) else {
            // 全体の係数が0になるので、評価結果も0になる
            continue;
        };
        let mut is_zero = false;
        let mut mahler_vec = Vec::new();
        for (&e, r) in exp.as_slice().iter().zip(sample.iter()) {
            if e == 0 {
                continue;
            }
            let Some(mahler) = get_signed_mahler(r) else {
                // このときは0が入っているので、全体として0になる
                is_zero = true;
                break;
            };
            mahler_vec.push(pow_signed_mahler(&mahler, e));
        }
        if is_zero {
            continue;
        }
        let term_mahler = mahler_vec
            .into_iter()
            .fold(coeff_mahler, |acc, m| mul_signed_mahler(&acc, &m));
        term_mahlers.push(term_mahler);
    }

    if let Some((first, rest)) = term_mahlers.split_first() {
        let first = signed_to_unsigned_mahler(first);
        let rest = rest
            .iter()
            .map(signed_to_unsigned_mahler)
            .collect::<Vec<_>>();
        // ここでfirstとrestを足し合わせる
        let result = rest
            .into_iter()
            .fold(first, |acc, m| add_unsigned_mahler(&acc, &m));
        if result.lower > BigRational::zero() {
            MahlerResult::Positive
        } else if result.upper < BigRational::zero() {
            MahlerResult::Negative
        } else {
            // 幅を計算する
            let width = &result.upper - &result.lower;
            // Landauの不等式による上界と比較する。width * 2^log_landau_bound < 1ならば符号が確定する
            let bound =
                BigRational::from_integer(BigInt::from(2).pow(result.log_landau_bound as u32));
            if width * bound < BigRational::from_integer(BigInt::from(1)) {
                // zeroになる
                MahlerResult::Zero
            } else {
                MahlerResult::Uncertain
            }
        }
    } else {
        // 全ての項が0だった場合は0になる
        MahlerResult::Zero
    }
}

impl MahlerResult {
    // MahlerResultがRelOpの等号を満たすかどうかを判定する関数
    pub fn satisfies(&self, rel_op: RelOp) -> Option<bool> {
        match self {
            MahlerResult::Uncertain => None,
            MahlerResult::Zero => Some(matches!(rel_op, RelOp::Eq | RelOp::Leq | RelOp::Geq)),
            MahlerResult::Positive => Some(matches!(rel_op, RelOp::Neq | RelOp::Gt | RelOp::Geq)),
            MahlerResult::Negative => Some(matches!(rel_op, RelOp::Neq | RelOp::Lt | RelOp::Leq)),
        }
    }
}
