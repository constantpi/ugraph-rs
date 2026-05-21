use num::BigInt;
use num_traits::{One, Zero};
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

#[cfg(test)]
mod tests {
    use std::os::unix::raw::gid_t;

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
}
