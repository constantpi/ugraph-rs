use num::BigInt;
use num_traits::One;
use vec1::Vec1;

use super::{BigIntPoly, PrimeModPoly, extended_gcd};

fn primepoly_to_bigintpoly(poly: &PrimeModPoly) -> BigIntPoly {
    let (rest, last) = poly.get_terms().clone().split_off_last();
    let rest = rest
        .iter()
        .cloned()
        .map(|c| c.to_bigint())
        .collect::<Vec<_>>();
    let last = last.to_bigint();
    let coeffs = Vec1::from_vec_push(rest, last);
    BigIntPoly::new(coeffs)
}

/// Henselの補題を用いて、f ≡ g * h (mod p) かつ gcd(g, h) = 1 を満たすときに、f ≡ G * H (mod p^k) かつ G ≡ g (mod p) かつ H ≡ h (mod p) を満たすG, Hを求める関数。
/// ただしhはmonicで、fの最高次の係数はgやhのprimeで割り切れないとする
fn lifting(
    f: &BigIntPoly,
    g: &PrimeModPoly,
    h: &PrimeModPoly,
    k: usize,
) -> (BigIntPoly, BigIntPoly) {
    let mut cnt = 0;
    let mut m = BigInt::from(g.get_terms().last().get_prime());
    let (s, t) = extended_gcd(g.clone(), h.clone());
    let mut s = primepoly_to_bigintpoly(&s);
    let mut t = primepoly_to_bigintpoly(&t);
    let mut g = primepoly_to_bigintpoly(g);
    let mut h = primepoly_to_bigintpoly(h);
    loop {
        if cnt == k {
            break (g, h);
        }
        cnt += 1;
        m = m.pow(2);

        let e = f.clone() - g.clone() * h.clone();
        let se = s.clone() * e.clone();
        let te = t.clone() * e.clone();
        let q = se.clone() / h.clone();
        let r = se - q.clone() * h.clone();
        let qg = q.clone() * g.clone();
        g = (g + te + qg).mod_integer(&m);
        h = (h + r).mod_integer(&m);

        let b = BigIntPoly::one() - s.clone() * g.clone() - t.clone() * h.clone();
        let sb = s.clone() * b.clone();
        let c = sb.clone() / h.clone();
        let d = sb - c.clone() * h.clone();
        let tb = t.clone() * b.clone();
        let cg = c.clone() * g.clone();
        s = (s + d).mod_integer(&m);
        t = (t + tb + cg).mod_integer(&m);
    }
}

fn lifting_vector(f: &BigIntPoly, g_vec: &Vec1<PrimeModPoly>) -> Vec1<BigIntPoly> {
    let p = g_vec.first().get_terms().last().get_prime();
    let len = g_vec.len_nonzero().get();
    if len == 1 {
        Vec1::new(f.clone())
    } else {
        // 長さが2以上のときは、g_vecを半分に分割して、それぞれをliftingしてから、さらにそれらをliftingする
        let mid = len / 2;
        let (first_half, second_half) = g_vec.split_at(mid);
        let g = first_half
            .iter()
            .fold(PrimeModPoly::one(p), |acc, g| acc * g.clone());
        let h = second_half
            .iter()
            .fold(PrimeModPoly::one(p), |acc, g| acc * g.clone());
        let (g_lifting, h_lifting) = lifting(f, &g, &h, 2);
        let mut g_vec_lifting = lifting_vector(
            &g_lifting,
            &Vec1::try_from_vec(first_half.to_vec()).unwrap(),
        );
        let h_vec_lifting = lifting_vector(
            &h_lifting,
            &Vec1::try_from_vec(second_half.to_vec()).unwrap(),
        );
        g_vec_lifting.extend(h_vec_lifting);
        g_vec_lifting
    }
}

#[cfg(test)]
mod tests {
    use super::super::PrimeField;
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_lifting() {
        let f = vec1![3.into(), 1.into(), 2.into(), 3.into()];
        let g = vec1![
            PrimeField::new(4, 5),
            PrimeField::new(1, 5),
            PrimeField::new(3, 5)
        ];
        let h = vec1![PrimeField::new(2, 5), PrimeField::new(1, 5)];
        let f = BigIntPoly::new(f);
        let g = PrimeModPoly::new(g, 5);
        let h = PrimeModPoly::new(h, 5);
        let (g_lifting, h_lifting) = lifting(&f, &g, &h, 2);
        let g_ans = BigIntPoly::new(vec1![209.into(), 301.into(), 3.into()]);
        let h_ans = BigIntPoly::new(vec1![317.into(), 1.into()]);
        assert_eq!(g_lifting, g_ans);
        assert_eq!(h_lifting, h_ans);
    }

    #[test]
    fn test_lifting_vector() {
        let f = vec1![1.into(), 3.into(), 3.into(), 2.into()];
        let f = BigIntPoly::new(f);
        let g_vec = vec1![
            PrimeModPoly::new(vec1![PrimeField::new(3, 7), PrimeField::new(2, 7)], 7),
            PrimeModPoly::new(vec1![PrimeField::new(3, 7), PrimeField::new(1, 7)], 7),
            PrimeModPoly::new(vec1![PrimeField::new(4, 7), PrimeField::new(1, 7)], 7),
        ];
        let g_vec_lifting = lifting_vector(&f, &g_vec);
        // g_vec_liftingの中身をすべてかけ合わせる
        let product_lifting = g_vec_lifting
            .iter()
            .fold(BigIntPoly::one(), |acc, g| acc * g.clone());
        let diff = f - product_lifting;
        let diff_mod = diff.mod_integer(&BigInt::from(2401));
        assert!(diff_mod.is_zero());
    }
}
