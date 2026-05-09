use crate::cad::matrix::{Matrix, generate_neg_list};
use itertools::{Itertools, iproduct};
use num::BigRational;

use crate::polynomial::{Exponent, Polynomial};

/// Polynomialから、最高次の項の次数と係数を取り出す関数。
fn collect_terms(poly: &Polynomial) -> Option<Vec<(u32, Exponent, BigRational)>> {
    poly.raw_iter()
        .map(|(exp, coeff)| {
            exp.split_last()
                .map(|(last, rest)| (last, rest, coeff.clone()))
        })
        .collect()
}

/// Polynomialから末尾の変数についての係数を取り出す関数
fn collect_last_variable_coefficients(poly: &Polynomial) -> Option<Vec<Polynomial>> {
    let terms = collect_terms(poly)?;
    let max_degree = terms.iter().map(|(last, _, _)| *last).max()?;
    let mut coeffs = vec![Polynomial::zero(); (max_degree + 1) as usize];
    for (last, rest, coeff) in terms {
        let index = last as usize;
        coeffs[index].add_term(rest, coeff);
    }

    Some(coeffs)
}

/// 係数の列から微分を計算する関数
fn differentiate_coefficients(coeffs: &[Polynomial]) -> Vec<Polynomial> {
    coeffs
        .iter()
        .enumerate()
        .skip(1)
        .map(|(i, coeff)| coeff.mul_rational(BigRational::from_integer(i.into())))
        .collect()
}

/// 多項式が末尾の変数について定数であるかどうかを判定する関数。
fn is_constant_for_last_variable(poly: &Polynomial) -> bool {
    if let Some(terms) = collect_terms(poly) {
        terms.iter().all(|(last, _, _)| *last == 0)
    } else {
        true
    }
}

fn psc(f: &[Polynomial], g: &[Polynomial], l: usize) -> Polynomial {
    let Some(f_degree) = f.len().checked_sub(1) else {
        return Polynomial::zero();
    };
    let Some(g_degree) = g.len().checked_sub(1) else {
        return Polynomial::zero();
    };
    if l >= f_degree || l >= g_degree {
        return Polynomial::zero();
    }
    let f_rows = g_degree - l;
    let g_rows = f_degree - l;
    let size = f_rows + g_rows;
    let mut matrix = Matrix::zero(size).unwrap();
    for i in 0..f_rows {
        for j in 0..=f_degree {
            if size <= i + j {
                break;
            }
            matrix.set(i, i + j, f[f_degree - j].clone());
        }
    }
    for i in 0..g_rows {
        for j in 0..=g_degree {
            if size <= i + j {
                break;
            }
            matrix.set(f_rows + i, i + j, g[g_degree - j].clone());
        }
    }

    let neg_list = generate_neg_list(size);
    matrix.determinant(&neg_list)
}

fn psc_list(f: &[Polynomial], g: &[Polynomial]) -> Vec<Polynomial> {
    let min_degree = f.len().min(g.len());
    (0..min_degree).map(|l| psc(f, g, l)).collect()
}

/// 多項式から射影を計算する関数
pub fn project_polynomial(polys: &[Polynomial]) -> Vec<Polynomial> {
    let (constants, non_constants): (Vec<Polynomial>, Vec<Polynomial>) = polys
        .iter()
        .cloned()
        .partition(is_constant_for_last_variable);
    let mut projected = Vec::new();
    for constant in constants {
        projected.push(constant.constant_term().unwrap());
    }

    let coeff_lists: Vec<Vec<Polynomial>> = non_constants
        .iter()
        .filter_map(collect_last_variable_coefficients)
        .collect();

    // proj1の計算
    for coeffs in coeff_lists.iter() {
        for coeff in coeffs {
            projected.push(coeff.clone());
        }
    }

    // proj2の計算
    for coeffs in coeff_lists.iter() {
        let len = coeffs.len();
        for i in 2..=len {
            let f = &coeffs[..i];
            let df = &differentiate_coefficients(f);
            let pscs = psc_list(f, df);
            projected.extend(pscs);
        }
    }

    // proj3の計算
    for (coeffs1, coeffs2) in coeff_lists.iter().tuple_combinations() {
        let len1 = coeffs1.len();
        let len2 = coeffs2.len();
        for (i, j) in iproduct!(2..=len1, 2..=len2) {
            let f = &coeffs1[..i];
            let g = &coeffs2[..j];
            let pscs = psc_list(f, g);
            projected.extend(pscs);
        }
    }

    // 重複を削除する
    let mut unique_projected = Vec::new();
    for p in projected {
        if !p.is_fully_constant() && !unique_projected.contains(&p) {
            unique_projected.push(p);
        }
    }

    unique_projected
}
