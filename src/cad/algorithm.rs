use color_eyre::Result;

use super::{
    Root, UnivariatePolynomial, calc_sample_points, find_unique_roots, is_possible_solution,
    lifting, polynomial_to_univariate, project_polynomial,
};
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
    let all_roots = find_unique_roots(&univariate_polynomials);
    for root in &all_roots {
        println!("Found root: {}", root);
    }
    let mut sample_points = calc_sample_points(&all_roots)
        .iter()
        .map(|r| vec![r.clone()])
        .collect::<Vec<_>>();
    println!("Sample points length: {}", sample_points.len());
    for sample in sample_points.iter() {
        println!("Trying sample point: {}", sample[0]);
    }
    // historyを逆順にたどりながら、sample_pointsをliftingしていく
    for polynomials in history.into_iter().rev() {
        sample_points = lifting(&polynomials, &sample_points)?;
        println!(
            "After lifting, sample points length: {}",
            sample_points.len()
        );
        for sample in &sample_points {
            print_sample_point(sample);
        }
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
            print_sample_point(&sample);
            if is_possible_solution(&polinomials, &sample)? {
                println!("Found possible solution");
                acc.push(sample);
            } else {
                println!("Sample point is not a solution");
            }
        }
        acc
    };
    for solution in &possible_solutions {
        print_sample_point(solution);
    }

    todo!()
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
