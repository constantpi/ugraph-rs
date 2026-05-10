use num::BigRational;
use vec1::Vec1;

use crate::cad::univariate::{
    UnivariatePolynomial, uni_poly_derivative, uni_poly_div, uni_poly_remainder,
};

/// オイラーの互除法
fn eulers_algorithm(poly: &UnivariatePolynomial) -> Vec1<UnivariatePolynomial> {
    let mut sequence = Vec1::new(poly.clone());
    let derivative = uni_poly_derivative(poly);
    sequence.push(derivative.clone());
    let (mut f, mut g) = (poly.clone(), derivative);
    loop {
        let r = uni_poly_remainder(&f, &g);
        if let Some(r) = r
            && !r.is_zero()
        {
            sequence.push(r.clone());
            f = g;
            g = r;
        } else {
            break sequence;
        }
    }
}

fn square_free_part(poly: &UnivariatePolynomial) -> UnivariatePolynomial {
    let deg = poly.degree();
    if deg <= 1 {
        poly.clone()
    } else {
        let sequence = eulers_algorithm(poly);
        let last = sequence.last();
        uni_poly_div(poly, last).unwrap().monic()
    }
}

pub fn sturm(poly: &UnivariatePolynomial) {
    let square_free = square_free_part(poly);
    let sequence = eulers_algorithm(&square_free);
}

#[cfg(test)]
mod tests {
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_square_free_part() {
        // (x+1)(x+1)(x-1) = x^3 + x^2 - x - 1 の平方因子を取り除くと x^2 - 1 になるはず
        let poly = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-1).into()),
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(1.into()),
            BigRational::from_integer(1.into())
        ]);
        let result = square_free_part(&poly);
        let expected = UnivariatePolynomial::new(vec1![
            BigRational::from_integer((-1).into()),
            BigRational::from_integer(0.into()),
            BigRational::from_integer(1.into())
        ]);
        assert_eq!(result, expected);
    }
}
