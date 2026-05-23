use num::{BigRational, Signed};
use num_traits::{One, Zero};

use super::{Range, Root, refine_range};
use crate::polynomial::Polynomial;

/// 符号の確定した区間を表す構造体
struct SignedInterval {
    is_positive: bool,
    abs_lower: BigRational,
    abs_upper: BigRational,
}

/// 符号の確定した区間または厳密値を表す構造体
enum SignedRange {
    Exact(BigRational),
    Interval(SignedInterval),
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
fn get_signed_range(root: &Root) -> SignedRange {
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
                let (new_low, new_high) = refine_range(poly, &lower, &upper);
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

fn pow_range(base: &SignedRange, exp: u32) -> SignedRange {
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

fn mul_ranges(range1: &SignedRange, range2: &SignedRange) -> SignedRange {
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

pub fn evaluate_polynomial_at_signed_range(
    poly: &Polynomial,
    sample: &[Root],
) -> (BigRational, BigRational) {
    let sample_ranges = sample.iter().map(get_signed_range).collect::<Vec<_>>();
    let mut lower = BigRational::zero();
    let mut upper = BigRational::zero();
    for (exp, coeff) in poly.raw_iter() {
        let product = exp
            .as_slice()
            .iter()
            .zip(sample_ranges.iter())
            .map(|(&e, r)| pow_range(r, e))
            .fold(SignedRange::Exact(BigRational::one()), |acc, r| {
                mul_ranges(&acc, &r)
            });
        let term_value = mul_constant_range(&product, coeff);
        match term_value {
            SignedRange::Exact(r) => {
                upper += r.clone();
                lower += r;
            }
            SignedRange::Interval(interval) => {
                if interval.is_positive {
                    upper += interval.abs_upper;
                    lower += interval.abs_lower;
                } else {
                    upper -= interval.abs_lower;
                    lower -= interval.abs_upper;
                }
            }
        }
    }
    (lower, upper)
}
