mod coordinate;
mod graph;
mod groebner;
mod polynomial;

use crate::coordinate::graph_to_polynomials;
use crate::graph::{generate_graph, simplify_graph};

use color_eyre::Result;

fn main() -> Result<()> {
    println!("Hello, world!");

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
