use num::BigRational;

use crate::graph::Graph;
use crate::groebner::groebner_basis;
use crate::polynomial::{Exponent, Polynomial};

/// xyを表す項を生成する関数
fn xy_term(n: usize, x: usize, y: usize) -> Exponent {
    let mut exponents = vec![0; n * 2];
    exponents[x] += 1; // x座標の変数の指数を1増やす
    exponents[y] += 1; // y座標の変数の指数を1増やす
    Exponent::new(exponents)
}

fn x_term(n: usize, x: usize) -> Exponent {
    let mut exponents = vec![0; n * 2];
    exponents[x] += 1; // x座標の変数の指数を1増やす
    Exponent::new(exponents)
}

fn constant_term(n: usize) -> Exponent {
    Exponent::new(vec![0; n * 2])
}

/// グラフを受け取って等長制約を表す多項式を生成する
pub fn graph_to_polynomials(graph: &Graph) -> Vec<Polynomial> {
    let mut polynomials = Vec::new();
    let n = graph.get_node_num();
    let one = BigRational::from_integer(1.into());
    let minus_one = BigRational::from_integer((-1).into());
    let minus_two = BigRational::from_integer((-2).into());
    for node in graph.iter_nodes() {
        let Some(neighbors) = graph.neighbors(node) else {
            continue;
        };
        for &neighbor in neighbors {
            if node > neighbor {
                continue; // 同じ辺を二回処理しないようにする
            }
            let mut p = Polynomial::zero();
            // (x_node - x_neighbor)^2 + (y_node - y_neighbor)^2 - 1 = 0 を表す多項式を生成
            let x_node = node.get_id() * 2; // x座標の変数ID
            let y_node = node.get_id() * 2 + 1; // y座標の変数ID
            let x_neighbor = neighbor.get_id() * 2; // x座標の変数ID
            let y_neighbor = neighbor.get_id() * 2 + 1; // y座標の変数ID
            p.add_term(xy_term(n, x_node, x_node), one.clone()); // x_node^2
            p.add_term(xy_term(n, y_node, y_node), one.clone()); // y_node^2
            p.add_term(xy_term(n, x_neighbor, x_neighbor), one.clone()); // x_neighbor^2
            p.add_term(xy_term(n, y_neighbor, y_neighbor), one.clone()); // y_neighbor^2
            p.add_term(xy_term(n, x_node, x_neighbor), minus_two.clone()); // -2 * x_node * x_neighbor
            p.add_term(xy_term(n, y_node, y_neighbor), minus_two.clone()); // -2 * y_node * y_neighbor
            p.add_term(constant_term(n), minus_one.clone()); // -1
            polynomials.push(p);
        }
    }
    // 最初の頂点を原点に、それに隣接する頂点をx軸上に固定するための多項式を追加
    if let Some(first_node) = graph.iter_nodes().next() {
        let mut px = Polynomial::zero();
        let x_first = first_node.get_id() * 2; // x座標の変数ID
        px.add_term(x_term(n, x_first), one.clone()); // x_first = 0
        let mut py = Polynomial::zero();
        let y_first = first_node.get_id() * 2 + 1; // y座標の変数ID
        py.add_term(x_term(n, y_first), one.clone()); // y_first = 0
        polynomials.push(px);
        polynomials.push(py);
        if let Some(neighbors) = graph.neighbors(first_node)
            && let Some(neighbor) = neighbors.iter().min()
        {
            let mut p = Polynomial::zero();
            let x_neighbor = neighbor.get_id() * 2; // x座標の変数ID
            p.add_term(x_term(n, x_neighbor), one.clone()); // x_neighbor = 1
            polynomials.push(p);
            let y_neighbor = neighbor.get_id() * 2 + 1; // y座標の変数ID
            let mut p = Polynomial::zero();
            p.add_term(x_term(n, y_neighbor), one.clone()); // y_neighbor = 0
            p.add_term(constant_term(n), minus_one.clone()); // -1
            polynomials.push(p);
        }
    }
    groebner_basis(&polynomials)
    // polynomials
}
