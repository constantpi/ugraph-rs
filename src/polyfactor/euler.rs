use vec1::Vec1;

use super::{PrimeField, PrimeModPoly};
/// a, bが与えられたときにax+by=1を満たすx, yを求める関数
/// 足し算、引き算、掛け算、割り算が定義された環に対して定義できる
pub fn extended_gcd(a: PrimeModPoly, b: PrimeModPoly) -> (PrimeModPoly, PrimeModPoly) {
    let prime = a.get_prime();
    if b.is_zero() {
        let coeff = a
            .to_constant()
            .unwrap_or_else(|| panic!("a is not a constant polynomial: {}", a));
        let inv = PrimeField::new(1, prime) / coeff;
        (
            PrimeModPoly::new(Vec1::new(inv), prime),
            PrimeModPoly::zero(prime),
        )
    } else {
        let (x1, y1) = extended_gcd(b.clone(), a.clone() % b.clone());
        let x = y1.clone();
        let y = x1 - (a / b) * y1;
        (x, y)
    }
}
