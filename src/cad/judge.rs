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

/// SignRelationからabsを取り出す関数
fn get_cross_zero_abs(sign_relation: SignRelation) -> Option<BigRational> {
    match sign_relation {
        SignRelation::CrossZero(abs) => Some(abs),
        _ => None,
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
            let Some(abs) = get_cross_zero_abs(evaluate_polynomial_at_signed_range(poly, &sample))
            else {
                return Solution::NoSolution; // 0をまたいでいない場合は解ではない
            };
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

    let ratio = BigRational::from_integer(7.into()) / BigRational::from_integer(10.into());
    for poly in poly_list {
        let Some(abs) = get_cross_zero_abs(evaluate_polynomial_at_signed_range(poly, &sample))
        else {
            return Solution::NoSolution; // 0をまたいでいない場合は解ではない
        };
        let mut refined_sample = sample.clone();
        for i in 0..10 {
            for root in refined_sample.iter_mut() {
                if let Some(refined) = refine_root(root) {
                    *root = refined;
                }
            }
            let Some(new_abs) =
                get_cross_zero_abs(evaluate_polynomial_at_signed_range(poly, &refined_sample))
            else {
                return Solution::NoSolution; // 0をまたいでいない場合は解ではない
            };
            if new_abs > abs.clone() * (ratio.pow(i)) {
                return Solution::NoSolution; // 誤差が十分に縮まっていない場合は解ではない
            }
        }
    }
    Solution::Exist(sample)
}
