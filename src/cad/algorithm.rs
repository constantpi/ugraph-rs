use color_eyre::Result;

use super::{Root, find_all_roots, polynomial_to_univariate, project_polynomial};
use crate::polynomial::Polynomial;

pub enum Solution {
    NoSolution,
    Exist(Vec<Root>),
}

pub fn find_solution(polinomials: &[Polynomial]) -> Result<Solution> {
    // 変数の数をまずは取得
    let Some(num_vars) = polinomials
        .iter()
        .find_map(|p| p.raw_iter().next().map(|(exp, _)| exp.len()))
    else {
        return Ok(Solution::Exist(vec![])); // どの変数も出てこない場合は常に解が存在する
    };
    let mut current_num_vars = num_vars;
    let mut current_polynomials = polinomials.to_vec();
    let mut history = Vec::new();
    while current_num_vars > 1 {
        history.push(current_polynomials.clone());
        current_num_vars -= 1;
        current_polynomials = project_polynomial(&current_polynomials);
    }
    let univariate_polynomials = current_polynomials
        .iter()
        .map(polynomial_to_univariate)
        .collect::<Result<Vec<_>>>()?;
    let all_roots = univariate_polynomials
        .iter()
        .flat_map(find_all_roots)
        .collect::<Vec<_>>();
    for root in all_roots {
        println!("Checking root: {}", root);
    }

    todo!()
}
