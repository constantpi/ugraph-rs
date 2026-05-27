use color_eyre::Result;
use num_traits::Zero;

use super::{MahlerResult, evaluate_polynomial_by_mahler};
use crate::cad::{Root, Solution, evaluate_polynomial_at_constants, refine_root};
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
pub fn is_possible_solution_by_resultant(
    poly_list: &[Polynomial],
    sample: &[Root],
) -> Result<bool> {
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
pub fn is_solution_by_interval(poly_list: &[Polynomial], sample: &[Root]) -> Solution {
    let mut sample = sample.to_vec();

    for poly in poly_list {
        loop {
            let mahler_result = evaluate_polynomial_by_mahler(poly, &sample);
            match mahler_result {
                MahlerResult::Positive | MahlerResult::Negative => return Solution::NoSolution, // 符号が確定している場合は解ではない
                MahlerResult::Zero => break, // 0であることが確定している場合はこれ以上精度を高めない
                MahlerResult::Uncertain => {
                    // 符号が確定していない場合は、サンプル点をさらに精密化する
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
