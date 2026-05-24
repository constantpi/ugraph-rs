use color_eyre::Result;
use num::BigRational;
use num_traits::Zero;

use crate::cad::{
    Root, SignRelation, Solution, evaluate_polynomial_at_constants,
    evaluate_polynomial_at_signed_range, refine_root,
};
use crate::polynomial::Polynomial;

/// サンプル点が多項式の解である可能性があるかどうかを判定する関数
fn is_possible_root(poly: &Polynomial, sample: &[Root]) -> Result<bool> {
    // 多項式にサンプル点を代入して、その値が0であるかどうかを判定する
    let values = sample
        .iter()
        .map(|root| root.get_poly())
        .cloned()
        .collect::<Vec<_>>();
    let resultant = evaluate_polynomial_at_constants(poly, &values).ok_or_else(|| {
        color_eyre::eyre::eyre!("Failed to evaluate polynomial at the given sample points")
    })?;
    Ok(resultant.is_zero())
}

/// サンプル点が与えられた多項式系の解であるかどうかを判定する関数
pub fn is_possible_solution(poly_list: &[Polynomial], sample: &[Root]) -> Result<bool> {
    let mut poly_iter = poly_list.iter();
    loop {
        let Some(poly) = poly_iter.next() else {
            break Ok(true);
        };
        if !is_possible_root(poly, sample)? {
            break Ok(false);
        }
    }
}

/// 区間演算によってサンプル点が多項式の解であるかどうかを判定する関数
pub fn is_possible_solution_interval(poly_list: &[Polynomial], sample: &[Root]) -> Solution {
    let mut sample = sample.to_vec();
    // 1/2^20まで誤差が縮まっているかを判定する
    let half = BigRational::from_integer(1.into()) / BigRational::from_integer(2.into());
    let threshold = half.pow(20);

    for poly in poly_list {
        loop {
            match evaluate_polynomial_at_signed_range(poly, &sample) {
                SignRelation::Negative | SignRelation::Positive => return Solution::NoSolution, // 0をまたいでいない場合は解ではない
                SignRelation::CrossZero(abs) => {
                    // 0をまたいでいる場合は精度を高める
                    if abs < threshold {
                        break; // 十分に精度が高まっている場合はこれ以上精度を高めない
                    }
                    for root in sample.iter_mut() {
                        if let Some(refined) = refine_root(root) {
                            *root = refined;
                        }
                    }
                }
            }
        }
    }
    Solution::Exist(sample)
}
