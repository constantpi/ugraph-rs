use color_eyre::Result;

use super::{
    Root, calc_sample_points, find_unique_roots, is_solution_by_interval, lifting,
    polynomial_to_univariate, project_polynomial,
};
use crate::parser::RelOp;
use crate::polynomial::Polynomial;

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Solution {
    NoSolution,
    Exist(Vec<Root>),
}

pub fn find_solution_equality(polys: &[Polynomial]) -> Result<Solution> {
    let ineqs = polys
        .iter()
        .map(|p| (p.clone(), RelOp::Eq))
        .collect::<Vec<_>>();
    find_solution(&ineqs)
}

pub fn find_solution(ineqs: &[(Polynomial, RelOp)]) -> Result<Solution> {
    let polinomials = ineqs.iter().map(|(p, _)| p.clone()).collect::<Vec<_>>();
    // 変数の数をまずは取得
    let Some(num_vars) = polinomials
        .iter()
        .find_map(|p| p.raw_iter().next().map(|(exp, _)| exp.len()))
    else {
        return Ok(Solution::Exist(vec![])); // どの変数も出てこない場合は常に解が存在する
    };
    let ineqs_by_num_vars = {
        let mut ineqs_by_num_vars = vec![Vec::new(); num_vars + 1];
        for (poly, rel_op) in ineqs {
            let poly_num_vars = poly.num_vars();
            if poly_num_vars > num_vars {
                return Err(color_eyre::eyre::eyre!(
                    "Unexpected error: Polynomial has more variables than expected"
                ));
            }
            ineqs_by_num_vars[poly_num_vars].push((poly.clone(), *rel_op));
        }
        ineqs_by_num_vars
    };
    // まず変数が0個のものは処理しておく
    if is_solution_by_interval(&ineqs_by_num_vars[0], &[]) == Solution::NoSolution {
        return Ok(Solution::NoSolution);
    }
    let mut current_num_vars = num_vars;
    let mut current_polynomials = polinomials.to_vec();
    let mut history = Vec::new();
    while current_num_vars > 1 {
        history.push(current_polynomials.clone());
        current_num_vars -= 1;
        current_polynomials = project_polynomial(&current_polynomials);
        println!(
            "Finished projection. Current number of variables: {}, number of projected polynomials: {}",
            current_num_vars,
            current_polynomials.len()
        );
    }
    let univariate_polynomials = current_polynomials
        .iter()
        .map(polynomial_to_univariate)
        .collect::<Result<Vec<_>>>()?;
    let all_roots = find_unique_roots(&univariate_polynomials);
    let mut sample_points = calc_sample_points(&all_roots)
        .iter()
        .map(|r| vec![r.clone()])
        .collect::<Vec<_>>();
    println!("Sample points after projection: {}", sample_points.len());
    // historyを逆順にたどりながら、sample_pointsをliftingしていく
    for polynomials in history.into_iter().rev() {
        // liftingの前に、sample_pointsがineqsを満たすかどうかを判定する
        sample_points.retain(|sample| {
                matches!(
                    is_solution_by_interval(&ineqs_by_num_vars[current_num_vars], sample),
                    Solution::Exist(_)
                )
            });
        println!(
            "sample points before lifting: {}, current_num_vars: {}",
            sample_points.len(),
            current_num_vars
        );
        sample_points = lifting(&polynomials, &sample_points)?;
        current_num_vars += 1;
        println!(
            "sample points after lifting: {}, current_num_vars: {}",
            sample_points.len(),
            current_num_vars
        );
    }
    sample_points.retain(|sample| {
            matches!(
                is_solution_by_interval(&ineqs_by_num_vars[current_num_vars], sample),
                Solution::Exist(_)
            )
        });
    if sample_points.iter().any(|sample| sample.len() != num_vars) {
        return Err(color_eyre::eyre::eyre!(
            "Unexpected error: Sample points have incorrect number of variables"
        ));
    }
    println!("Total sample points after lifting: {}", sample_points.len());
    let mut ans = None;
    for solution in &sample_points {
        if let Solution::Exist(refined_solution) = is_solution_by_interval(ineqs, solution) {
            println!("Sample point: {}", sample_point_to_string(solution));
            println!("This sample point is a solution.");
            ans = Some(refined_solution.clone());
        }
    }
    Ok(match ans {
        Some(solution) => Solution::Exist(solution),
        None => Solution::NoSolution,
    })
}

fn sample_point_to_string(sample: &[Root]) -> String {
    let mut s = String::new();
    s.push('(');
    for (i, root) in sample.iter().enumerate() {
        if i > 0 {
            s.push_str(", ");
        }
        s.push_str(&format!("{}", root));
    }
    s.push(')');
    s
}

impl std::fmt::Display for Solution {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Solution::NoSolution => write!(f, "No solution"),
            Solution::Exist(sample) => {
                write!(f, "Exist: (")?;
                for (i, root) in sample.iter().enumerate() {
                    if i > 0 {
                        write!(f, ", ")?;
                    }
                    write!(f, "{}", root)?;
                }
                write!(f, ")")
            }
        }
    }
}
