use color_eyre::Result;

use super::{
    Root, calc_sample_points, find_unique_roots, is_possible_solution_by_resultant,
    is_solution_by_interval, lifting, polynomial_to_univariate, project_polynomial,
};
use crate::polynomial::Polynomial;

#[derive(Debug, Clone, PartialEq, Eq)]
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
    println!(
        "Finished projection. Number of projected polynomials: {}",
        current_polynomials.len()
    );
    let univariate_polynomials = current_polynomials
        .iter()
        .map(polynomial_to_univariate)
        .collect::<Result<Vec<_>>>()?;
    let all_roots = find_unique_roots(&univariate_polynomials);
    let mut sample_points = calc_sample_points(&all_roots)
        .iter()
        .map(|r| vec![r.clone()])
        .collect::<Vec<_>>();
    // historyを逆順にたどりながら、sample_pointsをliftingしていく
    for polynomials in history.into_iter().rev() {
        println!("sample points before lifting: {}", sample_points.len());
        sample_points = lifting(&polynomials, &sample_points)?;
    }
    if sample_points.iter().any(|sample| sample.len() != num_vars) {
        return Err(color_eyre::eyre::eyre!(
            "Unexpected error: Sample points have incorrect number of variables"
        ));
    }
    println!("Total sample points after lifting: {}", sample_points.len());
    let possible_solutions = {
        let mut acc = Vec::new();
        for sample in sample_points {
            if is_possible_solution_by_resultant(polinomials, &sample)? {
                acc.push(sample);
            }
        }
        acc
    };
    println!(
        "Possible solutions after resultant check: {}",
        possible_solutions.len()
    );
    let mut ans = None;
    for solution in &possible_solutions {
        if let Solution::Exist(refined_solution) = is_solution_by_interval(polinomials, solution) {
            print_sample_point(solution);
            println!("This sample point is a possible solution.");
            ans = Some(refined_solution.clone());
        }
    }
    Ok(match ans {
        Some(solution) => Solution::Exist(solution),
        None => Solution::NoSolution,
    })
}

fn print_sample_point(sample: &[Root]) {
    print!("Sample point: (");
    for (i, root) in sample.iter().enumerate() {
        if i > 0 {
            print!(", ");
        }
        print!("{}", root);
    }
    println!(")");
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
