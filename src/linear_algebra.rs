use crate::polynomial::{Exponent, Polynomial};

use color_eyre::Result;
use num::BigRational;
use rand::Rng;
use rand::prelude::*;

/// 行列とベクトルからPolynomialを生成する関数
fn matrix_vector_to_polynomial(
    matrix: &Vec<Vec<BigRational>>,
    vector: &Vec<BigRational>,
) -> Result<Vec<Polynomial>> {
    if matrix.len() != vector.len() {
        return Err(color_eyre::eyre::eyre!(
            "The number of rows in the matrix must be equal to the length of the vector"
        ));
    }
    let mut polynomials = Vec::new();
    let mut len = None;
    for (row, b) in matrix.iter().zip(vector.iter()) {
        let mut p = Polynomial::zero();
        if let Some(l) = len
            && l != row.len()
        {
            return Err(color_eyre::eyre::eyre!(
                "All rows in the matrix must have the same length"
            ));
        }
        len = Some(row.len());
        for (i, a) in row.iter().enumerate() {
            let mut exponents = vec![0; row.len()];
            exponents[i] = 1; // x_iの指数を1にする
            p.add_term(Exponent::new(exponents), a.clone());
        }
        p.add_term(Exponent::new(vec![0; row.len()]), -b.clone()); // 定数項に-bを加える
        polynomials.push(p);
    }
    Ok(polynomials)
}

/// ランダムな行列とベクトルを生成する関数
fn generate_random_matrix_vector(n: usize, m: usize) -> (Vec<Vec<BigRational>>, Vec<BigRational>) {
    let mut rng = rand::rng();
    fn generate_random_int(rng: &mut impl rand::Rng) -> BigRational {
        let nums: Vec<i64> = (-10..=10).collect();
        let num = *nums.choose(rng).unwrap();
        BigRational::from_integer(num.into())
    }
    let matrix = (0..n)
        .map(|_| (0..m).map(|_| generate_random_int(&mut rng)).collect())
        .collect();
    let vector = (0..n).map(|_| generate_random_int(&mut rng)).collect();
    (matrix, vector)
}

pub fn generate_random_linear_polynomials(n: usize, m: usize) -> Result<Vec<Polynomial>> {
    let (matrix, vector) = generate_random_matrix_vector(n, m);
    matrix_vector_to_polynomial(&matrix, &vector)
}
