use itertools::Itertools;
use num::{BigInt, BigRational, Signed, Zero};

use super::determinant;
use crate::polynomial::{Exponent, Polynomial};

/// 多項式の指定された変数についての次数を計算する関数
fn degree_of_variable(poly: &Polynomial, var_index: usize) -> u32 {
    poly.raw_iter()
        .map(|(exp, _)| exp.get(var_index).unwrap_or(0))
        .max()
        .unwrap_or(0)
}

fn substitute_polynomial(poly: &Polynomial, subst_values: &[u32]) -> BigRational {
    poly.raw_iter()
        .map(|(exp, coeff)| {
            let mut term_value = coeff.clone();
            for (var_exp, value) in exp.iter().zip(subst_values.iter()) {
                term_value *= BigRational::from_integer(BigInt::from(*value).pow(var_exp));
            }
            term_value
        })
        .sum()
}

/// Ax = bとなるxを求める関数
fn solve_linear_system(a: &[Vec<BigRational>], b: &[BigRational]) -> Option<Vec<BigRational>> {
    let n = a.len();
    if a.iter().any(|row| row.len() != n) || b.len() != n {
        return None; // 正方行列でない場合やbのサイズが不適切な場合はNoneを返す
    }
    // ガウスの消去法を使って解を求める
    let mut mat = a.to_vec();
    let mut rhs = b.to_vec();
    for i in 0..n {
        // ピボット選択
        let mut pivot = i;
        for j in (i + 1)..n {
            if mat[j][i].abs() > mat[pivot][i].abs() {
                pivot = j;
            }
        }
        if mat[pivot][i].is_zero() {
            return None; // 解が存在しない場合はNoneを返す
        }
        if pivot != i {
            mat.swap(i, pivot);
            rhs.swap(i, pivot);
        }

        // 対策成分を1にする
        let pivot_value = mat[i][i].clone();
        for k in i..n {
            mat[i][k] = mat[i][k].clone() / pivot_value.clone();
        }
        rhs[i] = rhs[i].clone() / pivot_value;

        for j in 0..n {
            if i == j {
                continue;
            }
            let factor = mat[j][i].clone() / mat[i][i].clone();
            for k in i..n {
                mat[j][k] = mat[j][k].clone() - factor.clone() * mat[i][k].clone();
            }
            rhs[j] = rhs[j].clone() - factor.clone() * rhs[i].clone();
        }
    }

    Some(rhs)
}

/// Polynomialの行列の行列式を計算する関数。
/// ただし、変数があるとN!通りの計算が必要になるため最高次数を各変数について計算しておいて、実際に値を代入してから行列式を計算する。
/// これを複数回行うことで行列式の結果である多項式の係数を解く。
pub fn polynomial_matrix_determinant(matrix: Vec<Vec<Polynomial>>) -> Option<Polynomial> {
    let n = matrix.len();
    if matrix.iter().any(|row| row.len() != n) {
        return None; // 正方行列でない場合はNoneを返す
    }
    // 変数の数を計算する
    let num_variables = matrix
        .iter()
        .flat_map(|row| row.iter())
        .flat_map(|poly| poly.raw_iter().map(|(exp, _)| exp.len()))
        .max()
        .unwrap_or(0);
    let degrees = (0..num_variables)
        .map(|var| {
            // 各変数について行ごとの最高次数を計算する
            let degrees_row_sum = matrix
                .iter()
                .map(|row| {
                    row.iter()
                        .map(|poly| degree_of_variable(poly, var))
                        .max()
                        .unwrap_or(0)
                })
                .sum::<u32>();
            let degree_col_sum = (0..n)
                .map(|col| {
                    (0..n)
                        .map(|row| degree_of_variable(&matrix[row][col], var))
                        .max()
                        .unwrap_or(0)
                })
                .sum::<u32>();
            degree_col_sum.min(degrees_row_sum)
        })
        .collect::<Vec<_>>();
    // 係数行列とその値を計算
    let mut coefficient_matrix = Vec::new();
    let mut constant_terms = Vec::new();
    for subst_values in degrees
        .iter()
        .map(|&d| (0..=d).collect_vec())
        .multi_cartesian_product()
    {
        // subst_valuesは各変数に代入する値の組み合わせ
        let rational_matrix = matrix
            .iter()
            .map(|row| {
                row.iter()
                    .map(|poly| substitute_polynomial(poly, &subst_values))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        // rational_matrixの行列式を計算
        let det = determinant(&rational_matrix);
        let terms = degrees
            .iter()
            .map(|&d| (0..=d).collect_vec())
            .multi_cartesian_product()
            .map(|exps| {
                // expsは行列式の結果である多項式の項の次数の組み合わせ
                let term_value = exps
                    .iter()
                    .zip(subst_values.iter())
                    .map(|(&exp, &value)| BigRational::from_integer(BigInt::from(value).pow(exp)))
                    .product::<BigRational>();
                term_value
            })
            .collect::<Vec<_>>();
        coefficient_matrix.push(terms);
        constant_terms.push(det);
    }
    // coefficient_matrixとconstant_termsから多項式の係数を解く
    let coefficients = solve_linear_system(&coefficient_matrix, &constant_terms)?;
    // 係数と対応する項を組み合わせて多項式を構築する
    let mut result = Polynomial::zero();

    for (exps, coeff) in degrees
        .iter()
        .map(|&d| (0..=d).collect_vec())
        .multi_cartesian_product()
        .zip(coefficients)
    {
        let exps = Exponent::new(exps);
        result.add_term(exps, coeff);
    }
    Some(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_linear_system() {
        let a = vec![
            vec![
                BigRational::from_integer(2.into()),
                BigRational::from_integer(1.into()),
            ],
            vec![
                BigRational::from_integer(3.into()),
                BigRational::from_integer(2.into()),
            ],
        ];
        let b = vec![
            BigRational::from_integer(8.into()),
            BigRational::from_integer(3.into()),
        ];
        // 2x + y = 8
        // 3x + 2y = 3
        // 解はx = 13, y = -18
        let solution = solve_linear_system(&a, &b).unwrap();
        assert_eq!(solution[0], BigRational::from_integer(13.into()));
        assert_eq!(solution[1], BigRational::from_integer((-18).into()));
    }
}
