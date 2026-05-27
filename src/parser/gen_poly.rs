use color_eyre::Result;
use std::collections::HashMap;

use super::Expr;
use crate::polynomial::{Exponent, Polynomial};

/// まずはAstを見てそこから変数集合を特定する
fn extract_variables(expr: &Expr) -> Vec<String> {
    match expr {
        Expr::Num(_) => vec![],
        Expr::Var(name) => vec![name.clone()],
        Expr::Add(lhs, rhs) | Expr::Sub(lhs, rhs) | Expr::Mul(lhs, rhs) => {
            let mut vars = extract_variables(lhs);
            vars.extend(extract_variables(rhs));
            vars.sort();
            vars.dedup();
            vars
        }
        Expr::Neg(exp) | Expr::Pow(exp, _) => extract_variables(exp),
    }
}

fn gen_var(var_id: usize, var_num: usize) -> Polynomial {
    let mut p = Polynomial::zero();
    let mut exp_vec = vec![0; var_num];
    exp_vec[var_id] = 1;
    p.add_term(
        Exponent::new(exp_vec),
        num::BigRational::from_integer(1.into()),
    );
    p
}

fn gen_contant(c: i64, var_num: usize) -> Polynomial {
    let mut p = Polynomial::zero();
    p.add_term(
        Exponent::new(vec![0; var_num]),
        num::BigRational::from_integer(c.into()),
    );
    p
}

/// Astと変数へのマッピングから、Polynomialを生成する
fn ast_to_polynomial_using_map(
    expr: &Expr,
    var_map: &HashMap<String, usize>,
) -> Result<Polynomial> {
    let num_vars = var_map.len();
    let poly = match expr {
        Expr::Num(n) => gen_contant(*n, num_vars),
        Expr::Var(name) => {
            let var_id = var_map
                .get(name)
                .ok_or_else(|| color_eyre::eyre::eyre!("undefined variable: {}", name))?;
            gen_var(*var_id, num_vars)
        }
        Expr::Add(lhs, rhs) => {
            let lhs_poly = ast_to_polynomial_using_map(lhs, var_map)?;
            let rhs_poly = ast_to_polynomial_using_map(rhs, var_map)?;
            lhs_poly + rhs_poly
        }
        Expr::Sub(lhs, rhs) => {
            let lhs_poly = ast_to_polynomial_using_map(lhs, var_map)?;
            let rhs_poly = ast_to_polynomial_using_map(rhs, var_map)?;
            lhs_poly - rhs_poly
        }
        Expr::Mul(lhs, rhs) => {
            let lhs_poly = ast_to_polynomial_using_map(lhs, var_map)?;
            let rhs_poly = ast_to_polynomial_using_map(rhs, var_map)?;
            lhs_poly * rhs_poly
        }
        Expr::Neg(exp) => {
            let exp_poly = ast_to_polynomial_using_map(exp, var_map)?;
            -exp_poly
        }
        Expr::Pow(base, exp) => {
            let base_poly = ast_to_polynomial_using_map(base, var_map)?;
            if let Some(sub_exp) = (*exp).checked_sub(1) {
                (0..sub_exp).fold(base_poly.clone(), |acc, _| acc * base_poly.clone())
            } else {
                gen_contant(1, num_vars)
            }
        }
    };
    Ok(poly)
}

/// AstからPolynomialを生成する
pub fn ast_to_polynomial(expr: &[Expr]) -> Result<Vec<Polynomial>> {
    let vars = expr.iter().flat_map(extract_variables).collect::<Vec<_>>();
    // sortして重複を削除
    let vars = {
        let mut vars = vars;
        vars.sort();
        vars.dedup();
        vars
    };
    let var_map = vars
        .into_iter()
        .enumerate()
        .map(|(i, name)| (name, i))
        .collect();
    expr.iter()
        .map(|e| ast_to_polynomial_using_map(e, &var_map))
        .collect()
}
