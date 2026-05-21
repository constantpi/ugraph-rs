use itertools::Itertools;
use num::BigInt;
use vec1::Vec1;

use super::{
    PrimeField, PrimeModPoly, find_ok_prime, gcd, matrix_kernel, mod_poly_division,
    mod_poly_remainder,
};

/// 因数分解のためのBerlekampのアルゴリズム
pub fn berlekamp_factorization(coeffs: Vec1<BigInt>) -> Vec<PrimeModPoly> {
    let lt = coeffs.last().clone();
    let poly = find_ok_prime(coeffs);
    let p = poly.get_prime();
    let degree = poly.degree();
    let zero = PrimeField::zero(p);
    let lt: usize = (lt % p).try_into().unwrap();
    let lt = PrimeField::new(lt, p);

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

        let kernel_vectors = matrix_kernel(&matrix);
        // 不動点多項式
        let fixed_polys = kernel_vectors
            .into_iter()
            .map(|vec| {
                let coeffs = Vec1::try_from_vec(vec).unwrap();
                PrimeModPoly::new(coeffs, p)
            })
            .filter(|p| p.degree() > 0)
            .collect::<Vec<_>>();
        let mut ans: Vec<PrimeModPoly> = vec![];
        let mut queue = vec![poly];
        let candidates = fixed_polys
            .iter()
            .cartesian_product((0..p).map(|c| PrimeField::new(c, p)))
            .map(|(f, c)| f.add_const(&c))
            .collect::<Vec<_>>();
        loop {
            let Some(f) = queue.pop() else {
                if let Some(fist) = ans.first_mut() {
                    *fist = fist.mul_const(&lt);
                }
                break ans;
            };
            let f_degree = f.degree();
            if let Some(d) = candidates.iter().find_map(|g| {
                let d = gcd(&f, g);
                let d_degree = d.degree();
                if d_degree > 0 && d_degree < f_degree {
                    Some(d)
                } else {
                    None
                }
            }) {
                queue.push(mod_poly_division(&f, &d).unwrap());
                queue.push(d);
            } else {
                ans.push(f.monic());
            }
        }
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
        assert_eq!(factors.len(), 3);

        let coeffs = vec1![
            2.into(),
            2.into(),
            2.into(),
            4.into(),
            2.into(),
            4.into(),
            0.into(),
            2.into()
        ]; // (x^7 + 2x^5 + x^4 + 2x^3 + x^2 + x + 1) * 2

        let factors = berlekamp_factorization(coeffs);
        assert_eq!(factors.len(), 3);
    }
}
