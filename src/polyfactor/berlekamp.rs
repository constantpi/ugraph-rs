use num::BigInt;
use vec1::Vec1;

use super::{PrimeField, PrimeModPoly, find_ok_prime, mod_poly_remainder};

/// 因数分解のためのBerlekampのアルゴリズム
pub fn berlekamp_factorization(coeffs: Vec1<BigInt>) -> Vec<PrimeModPoly> {
    let poly = find_ok_prime(coeffs);
    let p = poly.get_prime();
    let degree = poly.degree();
    let zero = PrimeField::zero(p);
    if degree <= 1 {
        vec![poly]
    } else {
        let berlekamp_matrix = (0..degree)
            .map(|i| {
                // x^ip mod polyを計算する
                let coeffs = {
                    let zeros = vec![zero; i * p];
                    Vec1::from_vec_push(zeros, PrimeField::one(p))
                };
                let xip = PrimeModPoly::new(coeffs, p);
                let xip_mod = mod_poly_remainder(&xip, &poly)
                    .unwrap_or(PrimeModPoly::new(Vec1::new(zero), p));
                let terms = xip_mod.get_terms();
                (0..degree)
                    .map(|j| terms.get(j).unwrap_or(&zero))
                    .cloned()
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let matrix = (0..degree)
            .map(|i| {
                (0..degree)
                    .map(|j| {
                        if i == j {
                            berlekamp_matrix[j][i].clone() - PrimeField::one(p)
                        } else {
                            berlekamp_matrix[j][i].clone()
                        }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();

        for row in matrix.iter() {
            for coeff in row.iter() {
                print!("{} ", coeff);
            }
            println!();
        }
        todo!()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_berlekamp_factorization() {
        let coeffs = vec1![
            1.into(),
            1.into(),
            1.into(),
            2.into(),
            1.into(),
            2.into(),
            0.into(),
            1.into()
        ]; // x^7 + 2x^5 + x^4 + 2x^3 + x^2 + x + 1
        let factors = berlekamp_factorization(coeffs);
    }
}
