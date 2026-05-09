mod cad;
mod coordinate;
mod graph;
mod groebner;
mod polynomial;

use crate::cad::project_polynomial;
use crate::coordinate::graph_to_polynomials;
use crate::graph::{generate_graph, simplify_graph};

use color_eyre::Result;

fn main() -> Result<()> {
    println!("Hello, world!");

    let mut polynomial = polynomial::Polynomial::zero();
    polynomial.add_term(
        polynomial::Exponent::new(vec![2, 0, 0]),
        num::BigRational::from_integer((1).into()),
    ); // x^2
    polynomial.add_term(
        polynomial::Exponent::new(vec![0, 2, 0]),
        num::BigRational::from_integer((1).into()),
    ); // y^2
    polynomial.add_term(
        polynomial::Exponent::new(vec![0, 0, 2]),
        num::BigRational::from_integer((1).into()),
    ); // z^2
    polynomial.add_term(
        polynomial::Exponent::new(vec![1, 0, 0]),
        num::BigRational::from_integer((-4).into()),
    ); // -4x
    polynomial.add_term(
        polynomial::Exponent::new(vec![0, 1, 0]),
        num::BigRational::from_integer((-4).into()),
    ); // -4y
    polynomial.add_term(
        polynomial::Exponent::new(vec![0, 0, 1]),
        num::BigRational::from_integer((-4).into()),
    ); // -4z
    polynomial.add_term(
        polynomial::Exponent::new(vec![0, 0, 0]),
        num::BigRational::from_integer((11).into()),
    ); // +11
    println!("Original polynomial: {}", polynomial);
    println!("######################################################");
    let projected = project_polynomial(&[polynomial]);
    for (i, poly) in projected.iter().enumerate() {
        println!("Projected polynomial {}: {}", i + 1, poly);
    }
    println!("######################################################");
    let projected2 = project_polynomial(&projected);
    for (i, poly) in projected2.iter().enumerate() {
        println!("Projected polynomial {}: {}", i + 1, poly);
    }
    println!("######################################################");
    let f1 = {
        let mut p = polynomial::Polynomial::zero();
        p.add_term(
            polynomial::Exponent::new(vec![1, 1]),
            num::BigRational::from_integer(2431.into()),
        );
        p.add_term(
            polynomial::Exponent::new(vec![1, 0]),
            num::BigRational::from_integer((-2431).into()),
        );
        p.add_term(
            polynomial::Exponent::new(vec![0, 1]),
            num::BigRational::from_integer((-3301).into()),
        );
        p.add_term(
            polynomial::Exponent::new(vec![0, 0]),
            num::BigRational::from_integer(2685.into()),
        );
        p
    };
    let f2 = {
        let mut p = polynomial::Polynomial::zero();
        p.add_term(
            polynomial::Exponent::new(vec![0, 2]),
            num::BigRational::from_integer(1.into()),
        );
        p.add_term(
            polynomial::Exponent::new(vec![4, 0]),
            num::BigRational::from_integer(1.into()),
        );
        p.add_term(
            polynomial::Exponent::new(vec![1, 1]),
            num::BigRational::from_integer((-2).into()),
        );
        p
    };
    println!("f1: {}", f1);
    println!("f2: {}", f2);
    let projected = project_polynomial(&[f1, f2]);
    for (i, poly) in projected.iter().enumerate() {
        println!("Projected polynomial {}: {}", i + 1, poly);
    }
    println!("######################################################");

    let graph = generate_graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])?;
    for graph in simplify_graph(&graph)? {
        println!("Component:\n{}", graph);
        let polynomials = graph_to_polynomials(&graph);
        for (i, poly) in polynomials.iter().enumerate() {
            println!("f{}: {}", i + 1, poly);
        }
    }

    let graph = generate_graph(4, &[(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)])?;
    for graph in simplify_graph(&graph)? {
        println!("Component:\n{}", graph);
        let polynomials = graph_to_polynomials(&graph);
        for (i, poly) in polynomials.iter().enumerate() {
            println!("f{}: {}", i + 1, poly);
        }
    }

    let graph = generate_graph(
        7,
        &[
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (0, 5),
            (0, 6),
            (1, 2),
            (2, 3),
            (3, 4),
            (4, 5),
            (5, 6),
            (6, 1),
        ],
    )?;
    for graph in simplify_graph(&graph)? {
        println!("Component:\n{}", graph);
        let polynomials = graph_to_polynomials(&graph);
        for (i, poly) in polynomials.iter().enumerate() {
            println!("f{}: {}", i + 1, poly);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diamond() -> Result<()> {
        for _ in 0..10 {
            let graph = generate_graph(4, &[(0, 1), (1, 2), (2, 3), (3, 0), (0, 2)])?;
            let polynomials = graph_to_polynomials(&graph);
            assert_eq!(polynomials.len(), 11);
        }
        Ok(())
    }
}
