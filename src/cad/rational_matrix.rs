use num::{BigRational, One, Signed, Zero};

use crate::polynomial::Polynomial;

/// Polynomialの行列をBigRationalの行列に変換する関数
fn polynomial_matrix_to_rational_matrix(
    matrix: Vec<Vec<Polynomial>>,
) -> Option<Vec<Vec<BigRational>>> {
    matrix
        .iter()
        .map(|row| {
            row.iter()
                .map(|poly| poly.as_constant())
                .collect::<Option<Vec<_>>>()
        })
        .collect::<Option<Vec<_>>>()
}

/// BigRationalの行列の行列式を計算する関数
fn determinant(matrix: &[Vec<BigRational>]) -> BigRational {
    let n = matrix.len();
    if n == 0 {
        return BigRational::one();
    } else if n == 1 {
        return matrix[0][0].clone();
    }
    // n!通りの計算をすると非常に非効率なので、ガウスの消去法を使う
    let mut mat = matrix.to_vec();
    let mut det = BigRational::one();
    for i in 0..n {
        // ピボット選択
        let mut pivot = i;
        for j in (i + 1)..n {
            if mat[j][i].abs() > mat[pivot][i].abs() {
                pivot = j;
            }
        }
        if mat[pivot][i].is_zero() {
            return BigRational::from_integer(0.into());
        }
        if pivot != i {
            mat.swap(i, pivot);
            det = -det;
        }
        det *= mat[i][i].clone();
        for j in (i + 1)..n {
            let factor = mat[j][i].clone() / mat[i][i].clone();
            // for k in i..n {
            //     mat[j][k] = mat[j][k].clone() - factor.clone() * mat[i][k].clone();
            // }
            let row_i_tail: Vec<BigRational> = mat[i].iter().skip(i).cloned().collect();
            for (a, b) in mat[j].iter_mut().skip(i).zip(row_i_tail.into_iter()) {
                *a = (*a).clone() - factor.clone() * b;
            }
        }
    }
    det
}

/// Polynomialの行列の行列式を計算する関数
pub fn rational_matrix_determinant(matrix: Vec<Vec<Polynomial>>) -> Option<BigRational> {
    let rational_matrix = polynomial_matrix_to_rational_matrix(matrix)?;
    Some(determinant(&rational_matrix))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determinant() {
        let mat = vec![
            vec![
                BigRational::from_integer(1.into()),
                BigRational::from_integer(2.into()),
            ],
            vec![
                BigRational::from_integer(3.into()),
                BigRational::from_integer(4.into()),
            ],
        ];
        let det = determinant(&mat);
        assert_eq!(det, BigRational::from_integer((-2).into()));
    }
}
