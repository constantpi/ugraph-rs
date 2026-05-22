use super::PrimeField;

type Matrix = Vec<Vec<PrimeField>>;

/// Matrixを行基本変形によって階段行列に変換する関数。
/// 正方行列を想定している
fn to_row_echelon_form(matrix: &mut Matrix) {
    let n = matrix.len();
    for i in 0..n {
        // i列目のピボットを探す
        let Some(pivot_row) = (i..n).find(|&j| !matrix[j][i].is_zero()) else {
            // ピボットが見つからない場合は次の列に進む
            continue;
        };
        // ピボット行をi行目と交換する
        matrix.swap(i, pivot_row);
        let pivot = matrix[i][i];
        // ピボット行をpivotで割る
        matrix[i] = matrix[i].iter().map(|x| *x / pivot).collect();
        // i行目以外の行からi列目の成分を消す
        for j in 0..n {
            if j != i {
                let factor = matrix[j][i];
                matrix[j] = matrix[j]
                    .iter()
                    .zip(matrix[i].iter())
                    .map(|(x, y)| *x - factor * *y)
                    .collect();
            }
        }
    }
}

/// 行列のkernelを求める関数
pub fn matrix_kernel(matrix: &Matrix) -> Vec<Vec<PrimeField>> {
    let mut mat = matrix.clone();
    to_row_echelon_form(&mut mat);
    let free_vars = mat
        .iter()
        .enumerate()
        .filter_map(|(i, row)| if row[i].is_zero() { Some(i) } else { None });
    free_vars
        .map(|free_var| {
            mat.iter()
                .enumerate()
                .map(|(i, row)| {
                    if i == free_var {
                        PrimeField::one(row[0].get_prime())
                    } else {
                        -row[free_var]
                    }
                })
                .collect::<Vec<_>>()
        })
        .collect()
}
