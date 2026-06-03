use super::evaluate_polynomial_by_mahler;
use crate::cad::{Root, Solution, refine_root};
use crate::parser::RelOp;
use crate::polynomial::Polynomial;

/// 区間演算によってサンプル点が多項式の解であるかどうかを判定する関数
pub fn is_solution_by_interval(ineqs: &[(Polynomial, RelOp)], sample: &[Root]) -> Solution {
    let mut sample = sample.to_vec();

    for (poly, rel_op) in ineqs {
        loop {
            let mahler_result = evaluate_polynomial_by_mahler(poly, &sample);
            match mahler_result.satisfies(*rel_op) {
                None => {
                    // 符号が確定していない場合は、サンプル点をさらに精密化する
                    for root in sample.iter_mut() {
                        if let Some(refined) = refine_root(root, true) {
                            *root = refined;
                        }
                    }
                }
                Some(true) => break, // 条件を満たしていることが確定している場合はこれ以上精度を高めない
                Some(false) => return Solution::NoSolution, // 条件を満たしていないことが確定している場合は解ではない
            }
        }
    }

    Solution::Exist(sample)
}
