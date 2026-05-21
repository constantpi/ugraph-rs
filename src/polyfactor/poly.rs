use color_eyre::eyre::Ok;
use num::{BigInt, Zero};
use vec1::Vec1;

use super::{PrimeField, PrimeIter, is_prime};

fn clean(coefficients: &mut Vec1<PrimeField>) {
    // popが成功すれば残りもVec1であるため続行する
    if coefficients.last().is_zero() && coefficients.pop().is_ok() {
        clean(coefficients);
    }
}

/// Mod pでの多項式を表す構造体
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PrimeModPoly {
    // 項の指数と係数のマップ
    terms: Vec1<PrimeField>,
    prime: usize,
}

impl PrimeModPoly {
    pub fn new(terms: Vec1<PrimeField>, prime: usize) -> Self {
        if !is_prime(prime) {
            panic!("prime must be a prime number");
        }
        for term in &terms {
            if term.get_prime() != prime {
                panic!("All terms must have the same prime");
            }
        }
        let terms = {
            let mut t = terms;
            clean(&mut t);
            t
        };
        Self { terms, prime }
    }

    pub fn degree(&self) -> usize {
        self.terms.len_nonzero().get() - 1
    }

    pub fn iter(&self) -> impl Iterator<Item = &PrimeField> {
        self.terms.iter()
    }

    pub fn is_zero(&self) -> bool {
        self.terms.iter().all(|coeff| coeff.is_zero())
    }

    pub fn lt(&self) -> &PrimeField {
        self.terms.last()
    }

    pub fn get_prime(&self) -> usize {
        self.prime
    }

    pub fn get_terms(&self) -> &Vec1<PrimeField> {
        &self.terms
    }

    pub fn add_const(&self, constant: &PrimeField) -> Self {
        let mut new_terms = self.terms.clone();
        let first = new_terms.first_mut();
        *first = *first + *constant;
        PrimeModPoly::new(new_terms, self.prime)
    }

    pub fn mul_const(&self, constant: &PrimeField) -> Self {
        let (rest, last) = self.terms.clone().split_off_last();
        let new_last = last * *constant;
        let new_rest = rest.iter().map(|c| *c * *constant).collect::<Vec<_>>();
        let new_terms = Vec1::from_vec_push(new_rest, new_last);
        PrimeModPoly::new(new_terms, self.prime)
    }

    pub fn monic(&self) -> Self {
        let leading_coeff = self.lt();
        if leading_coeff.is_zero() {
            self.clone()
        } else {
            let (rest, last) = self.terms.clone().split_off_last();
            let rest_inv = rest.iter().map(|c| *c / last).collect::<Vec<_>>();
            let term_inv = Vec1::from_vec_push(rest_inv, PrimeField::one(self.prime));
            PrimeModPoly::new(term_inv, self.prime)
        }
    }

    pub fn one(prime: usize) -> Self {
        PrimeModPoly::new(Vec1::new(PrimeField::one(prime)), prime)
    }

    pub fn zero(prime: usize) -> Self {
        PrimeModPoly::new(Vec1::new(PrimeField::zero(prime)), prime)
    }

    pub fn to_constant(&self) -> Option<PrimeField> {
        if self.degree() == 0 {
            Some(self.terms.first().clone())
        } else {
            None
        }
    }
}

pub fn mod_poly_derivative(poly: &PrimeModPoly) -> PrimeModPoly {
    let new_terms = poly
        .iter()
        .enumerate()
        .skip(1) // 定数項は微分すると0になるためスキップ
        .map(|(i, coeff)| {
            *coeff * PrimeField::new(i as usize, poly.prime) // i*x^(i-1)の係数はi*coeffになる
        })
        .collect::<Vec<_>>();
    if let Some(coeffs) = Vec1::try_from_vec(new_terms).ok() {
        PrimeModPoly::new(coeffs, poly.prime)
    } else {
        // 全ての項が0になった場合は0多項式を返す
        PrimeModPoly::new(Vec1::new(PrimeField::new(0, poly.prime)), poly.prime)
    }
}

/// 多項式同士の割り算を行って余りを求める
/// ただしcleanがされていることが前提
fn remainder(f: &Vec1<PrimeField>, g: &Vec1<PrimeField>) -> Option<Vec1<PrimeField>> {
    let f_degree = f.len_nonzero().get() - 1;
    let g_degree = g.len_nonzero().get() - 1;
    let lt_f = f.last();
    let lt_g = g.last();
    let p = lt_g.get_prime();
    if g_degree == 0 {
        if lt_g.is_zero() {
            None // 0で割ることはできないためNoneを返す
        } else {
            Some(Vec1::new(PrimeField::zero(p))) // 定数で割るとあまりは常に0になる
        }
    } else if f_degree < g_degree {
        Some(f.clone()) // fの方が次数が小さい場合、割り切れないためfが余りになる
    } else {
        let mut ans = f.clone();
        // まずansの最高次の項を消す
        let leading_coeff = *lt_f / *lt_g;
        ans[f_degree] = PrimeField::zero(p);
        for i in 0..g_degree {
            let g_coeff = g[i];
            ans[f_degree - g_degree + i] -= leading_coeff * g_coeff;
        }
        // 末尾の0を削除する
        clean(&mut ans);
        remainder(&ans, g)
    }
}

fn div(f: &Vec1<PrimeField>, g: &Vec1<PrimeField>) -> Option<Vec1<PrimeField>> {
    let f_degree = f.len_nonzero().get() - 1;
    let g_degree = g.len_nonzero().get() - 1;
    let lt_f = f.last();
    let lt_g = g.last();
    let p = lt_g.get_prime();
    if lt_g.is_zero() {
        None // 0で割ることはできないためNoneを返す
    } else {
        match f_degree.cmp(&g_degree) {
            std::cmp::Ordering::Less => Some(Vec1::new(PrimeField::zero(p))), // fの方が次数が小さい場合、商は常に0になる
            std::cmp::Ordering::Equal => Some(Vec1::new(*lt_f / *lt_g)), // 同じ次数の場合、商は定数になる
            std::cmp::Ordering::Greater => {
                let res = *lt_f / *lt_g;
                let mut f = f.clone();
                f.pop().unwrap(); // 最高次の項を削除。fの次数はg以上なので必ず成功する
                for i in 0..g_degree {
                    let g_coeff = g[i].clone();
                    f[f_degree - g_degree + i] -= res.clone() * g_coeff;
                }
                let mut ans = div(&f, g)?;
                ans.push(res);
                Some(ans)
            }
        }
    }
}

pub fn mod_poly_division(f: &PrimeModPoly, g: &PrimeModPoly) -> Option<PrimeModPoly> {
    if f.prime != g.prime {
        panic!("Cannot divide polynomials over different fields");
    }
    div(&f.terms, &g.terms).map(|terms| PrimeModPoly::new(terms, f.prime))
}

pub fn mod_poly_remainder(f: &PrimeModPoly, g: &PrimeModPoly) -> Option<PrimeModPoly> {
    if f.prime != g.prime {
        panic!("Cannot divide polynomials over different fields");
    }
    remainder(&f.terms, &g.terms).map(|terms| PrimeModPoly::new(terms, f.prime))
}

pub fn gcd(f: &PrimeModPoly, g: &PrimeModPoly) -> PrimeModPoly {
    if f.prime != g.prime {
        panic!("Cannot compute gcd of polynomials over different fields");
    }
    let rem = mod_poly_remainder(f, g);
    if let Some(r) = rem {
        if r.is_zero() { g.clone() } else { gcd(g, &r) }
    } else {
        // gが0多項式の場合はfが最大公約数になる
        f.clone()
    }
}

fn is_ok_prime(coeffs: Vec1<BigInt>, p: usize) -> Option<PrimeModPoly> {
    let (rest, last) = coeffs.split_off_last();
    let last: usize = (last % p).try_into().unwrap();
    if last == 0 {
        // 最高次の項の係数がpで割り切れる場合はNG
        None
    } else {
        let last = PrimeField::new(last, p);
        let rest = rest
            .into_iter()
            .map(|c| PrimeField::new((c % p).try_into().unwrap(), p))
            .collect::<Vec<_>>();
        let mod_coeffs = Vec1::from_vec_push(rest, last);
        let poly = PrimeModPoly::new(mod_coeffs, p);
        let der_poly = mod_poly_derivative(&poly);
        let g = gcd(&poly, &der_poly);
        if g.degree() == 0 { Some(poly) } else { None }
    }
}

pub fn find_ok_prime(coeffs: Vec1<BigInt>) -> PrimeModPoly {
    let mut prime_iter = PrimeIter::new();
    loop {
        let p = prime_iter.next().unwrap();
        if let Some(poly) = is_ok_prime(coeffs.clone(), p) {
            break poly;
        }
    }
}

impl std::ops::Add for PrimeModPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.prime != rhs.prime {
            panic!("Cannot add polynomials over different fields");
        }
        let (constant, non_constant) = self.terms.split_off_first();
        let (rhs_constant, rhs_non_constant) = rhs.terms.split_off_first();
        let mut sum = Vec1::new(constant + rhs_constant);
        let max_len = usize::max(non_constant.len(), rhs_non_constant.len());
        for i in 0..max_len {
            let coeff1 = non_constant
                .get(i)
                .cloned()
                .unwrap_or(PrimeField::zero(self.prime));
            let coeff2 = rhs_non_constant
                .get(i)
                .cloned()
                .unwrap_or(PrimeField::zero(self.prime));
            sum.push(coeff1 + coeff2);
        }
        PrimeModPoly::new(sum, self.prime)
    }
}

impl std::ops::Neg for PrimeModPoly {
    type Output = Self;

    fn neg(self) -> Self::Output {
        let (constant, non_constant) = self.terms.split_off_first();
        let mut neg_coeffs = Vec1::new(-constant);
        for coeff in non_constant {
            neg_coeffs.push(-coeff);
        }
        PrimeModPoly::new(neg_coeffs, self.prime)
    }
}

impl std::ops::Sub for PrimeModPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl std::ops::Mul for PrimeModPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.prime != rhs.prime {
            panic!("Cannot multiply polynomials over different fields");
        }
        let p = self.prime;
        let mut new_coeffs = Vec1::new(PrimeField::zero(p));
        new_coeffs.extend(vec![PrimeField::zero(p); self.degree() + rhs.degree()]);
        for (i, coeff1) in self.terms.iter().enumerate() {
            for (j, coeff2) in rhs.terms.iter().enumerate() {
                new_coeffs[i + j] += *coeff1 * *coeff2;
            }
        }
        PrimeModPoly::new(new_coeffs, p)
    }
}

impl std::ops::Div for PrimeModPoly {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        mod_poly_division(&self, &rhs).expect("Division by zero polynomial")
    }
}

impl std::ops::Rem for PrimeModPoly {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self::Output {
        mod_poly_remainder(&self, &rhs).expect("Division by zero polynomial")
    }
}

impl std::fmt::Display for PrimeModPoly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let degree = self.degree();
        for (i, coeff) in self.terms.iter().enumerate().rev() {
            let exp = if i == 0 {
                String::new()
            } else if i == 1 {
                "x".to_string()
            } else {
                format!("x^{}", i)
            };
            let coeff = if coeff.is_zero() {
                continue;
            } else if coeff == &PrimeField::one(self.prime) && i != 0 {
                String::new()
            } else {
                format!("{}", coeff)
            };
            if i == degree {
                write!(f, "{}{}", coeff, exp)?;
            } else {
                write!(f, " + {}{} ", coeff, exp)?;
            }
        }
        write!(f, " (mod {})", self.prime)
    }
}

#[cfg(test)]
mod tests {
    use vec1::vec1;

    use super::super::extended_gcd;
    use super::*;

    #[test]
    fn test_is_ok_prime() {
        let coeffs = vec1![1.into(), 0.into(), 0.into(), 5.into()]; // 5x^3 + 1
        let result = is_ok_prime(coeffs, 5);
        assert!(result.is_none());

        let coeffs = vec1![3.into(), 0.into(), 1.into()]; // x^2 + 3
        let result = is_ok_prime(coeffs.clone(), 3);
        assert!(result.is_none());

        let result = is_ok_prime(coeffs, 5);
        assert!(result.is_some());
    }

    #[test]
    fn test_euler_gcd() {
        let f = vec1![
            PrimeField::new(4, 5),
            PrimeField::new(1, 5),
            PrimeField::new(3, 5)
        ];
        let g = vec1![PrimeField::new(2, 5), PrimeField::new(1, 5)];
        let f = PrimeModPoly::new(f, 5);
        let g = PrimeModPoly::new(g, 5);
        let (x, y) = extended_gcd(f.clone(), g.clone());
        assert_eq!(f * x + g * y, PrimeModPoly::one(5));
    }
}
