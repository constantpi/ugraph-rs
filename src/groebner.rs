use num::BigRational;
use num::Zero;

use crate::parser::RelOp;
use crate::polynomial::Exponent;
use crate::polynomial::Polynomial;

// Polynomial構造体は係数が0の項を持たないことが保証されている。

/// 複数の多項式による割り算
fn multi_division(dividend: &Polynomial, divisors: &[Polynomial]) -> Polynomial {
    // ここに割り算のアルゴリズムを実装
    let mut remainder = Polynomial::zero();
    let mut p = dividend.clone();

    while let Some((exponent, coefficient)) = p.get_lt() {
        let exponent = exponent.clone();
        let coefficient = coefficient.clone();
        let mut divided = false;
        for divisor in divisors {
            if let Some((div_exponent, div_coefficient)) = divisor.get_lt()
                && let Some(sub) = exponent.sub(div_exponent)
            {
                // 割り算が可能な場合
                let factor = coefficient.clone() / div_coefficient;
                let term = divisor.mul_term(sub, factor);
                p = p - term;
                divided = true;
                break;
            }
        }
        if !divided {
            // どの割り算もできない場合、現在の項を剰余に追加
            remainder.add_term(exponent.clone(), coefficient.clone());
            p.add_term(exponent, -coefficient);
        }
    }

    remainder
}

fn s_polynomial(f: &Polynomial, g: &Polynomial) -> Option<Polynomial> {
    if let (Some((f_exp, f_coeff)), Some((g_exp, g_coeff))) = (f.get_lt(), g.get_lt()) {
        let f_term = f.mul_term(g_exp.saturated_sub(f_exp), g_coeff.clone());
        let g_term = g.mul_term(f_exp.saturated_sub(g_exp), f_coeff.clone());
        Some(f_term - g_term)
    } else {
        None
    }
}

pub fn groebner_basis(polynomials: &[Polynomial]) -> Vec<Polynomial> {
    let mut basis = polynomials.to_vec();
    let mut pairs = generate_pairs(&basis);
    loop {
        let Some((i, j)) = pairs.pop() else {
            // すべてのペアを処理し終えたら終了
            break basis;
        };
        if let Some(s) = s_polynomial(&basis[i], &basis[j]) {
            let r = multi_division(&s, &basis);
            if !r.is_zero() {
                basis.push(r);
                basis = simplify_groebner_basis(&basis);
                break groebner_basis(&basis);
            }
        }
    }
}
fn simplify_groebner_basis(polynomials: &[Polynomial]) -> Vec<Polynomial> {
    let mut polynomials = polynomials.to_vec();
    loop {
        let (changed, reduced) = remainder_reduce_basis(&polynomials);
        let len_before = polynomials.len();
        polynomials = remove_zero_polynomials(&reduced);
        polynomials = polynomials.into_iter().map(make_monic).collect();
        if !changed && polynomials.len() == len_before {
            break polynomials;
        }
    }
}

/// 優先度順のpairsを生成する
fn generate_pairs(basis: &[Polynomial]) -> Vec<(usize, usize)> {
    let mut pairs = Vec::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            let lcm_degree = if let (Some((f_exp, _)), Some((g_exp, _))) =
                (basis[i].get_lt(), basis[j].get_lt())
            {
                f_exp.lcm(g_exp).sum_degree()
            } else {
                0
            };
            pairs.push(((i, j), lcm_degree));
        }
    }

    // pairsをlcmの大きい順にソートする
    pairs.sort_by_key(|&(_, lcm_degree)| std::cmp::Reverse(lcm_degree));
    pairs.into_iter().map(|(pair, _)| pair).collect()
}

/// モニック化する
fn make_monic(poly: Polynomial) -> Polynomial {
    if let Some((exponent, coeff)) = poly.get_lt() {
        if coeff.is_zero() {
            poly.clone()
        } else {
            let inv_coeff = BigRational::from_integer(1.into()) / coeff.clone();
            poly.mul_term(Exponent::new(vec![0; exponent.len()]), inv_coeff)
        }
    } else {
        poly.clone()
    }
}
/// グレブナー基底から0多項式を削除する
fn remove_zero_polynomials(basis: &[Polynomial]) -> Vec<Polynomial> {
    basis.iter().filter(|p| !p.is_zero()).cloned().collect()
}

/// 基底の各要素を互いに割り算して、余りだけを残すことで基底をさらに簡略化する
fn remainder_reduce_basis(basis: &[Polynomial]) -> (bool, Vec<Polynomial>) {
    let mut reduced = basis.to_vec();
    let mut finished = false;
    let mut changed = false;
    while !finished {
        finished = true;
        for i in 0..reduced.len() {
            let r = multi_division(
                &reduced[i],
                &reduced
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != i)
                    .map(|(_, p)| p.clone())
                    .collect::<Vec<_>>(),
            );
            if r != reduced[i] {
                reduced[i] = r;
                finished = false;
                changed = true;
            }
        }
    }

    (changed, reduced)
}

/// 等式制約についてはグレブナー基底を求める
pub fn groebner_basis_for_equalities(
    polynomials: &[(Polynomial, RelOp)],
) -> Vec<(Polynomial, RelOp)> {
    let (equalities, inequalities): (Vec<_>, Vec<_>) = polynomials
        .iter()
        .cloned()
        .partition(|(_, op)| *op == RelOp::Eq);
    let basis = groebner_basis(&equalities.into_iter().map(|(p, _)| p).collect::<Vec<_>>());
    let mut result = basis
        .into_iter()
        .map(|p| (p, RelOp::Eq))
        .collect::<Vec<_>>();
    result.extend(inequalities.into_iter());
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::polynomial::{Exponent, Polynomial};
    use num::BigRational;

    #[test]
    fn test_groebner_basis() {
        let f1 = {
            let mut p = Polynomial::zero();
            p.add_term(
                Exponent::new(vec![1, 0, 0]),
                BigRational::from_integer(1.into()),
            );
            p.add_term(
                Exponent::new(vec![0, 1, 0]),
                BigRational::from_integer(1.into()),
            );
            p
        };
        let f2 = {
            let mut p = Polynomial::zero();
            p.add_term(
                Exponent::new(vec![1, 0, 0]),
                BigRational::from_integer(1.into()),
            );
            p.add_term(
                Exponent::new(vec![0, 0, 1]),
                BigRational::from_integer(1.into()),
            );
            p
        };
        let basis = groebner_basis(&[f1, f2]);

        assert_eq!(basis.len(), 2);
    }
}
