use color_eyre::Result;
use num::{BigRational, One, Zero};
use vec1::Vec1;

use crate::polynomial::Polynomial;

fn clean(coefficients: &mut Vec1<BigRational>) {
    // popが成功すれば残りもVec1であるため続行する
    if coefficients.last().is_zero() && coefficients.pop().is_ok() {
        clean(coefficients);
    }
}

/// 一変数多項式を表す構造体
/// i番目はx^iの係数を表す
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct UnivariatePolynomial(Vec1<BigRational>);

impl UnivariatePolynomial {
    pub fn new(coefficients: Vec1<BigRational>) -> Self {
        let mut ans = UnivariatePolynomial(coefficients);
        ans.clean();
        ans
    }

    pub fn degree(&self) -> usize {
        self.0.len_nonzero().get() - 1
    }

    pub fn iter(&self) -> impl Iterator<Item = &BigRational> {
        self.0.iter()
    }

    /// 末尾の0を削除する関数
    pub fn clean(&mut self) {
        clean(&mut self.0);
    }

    pub fn is_zero(&self) -> bool {
        self.0.iter().all(|coeff| coeff.is_zero())
    }

    pub fn mul_constant(&self, constant: &BigRational) -> Self {
        let (constant_coeff, non_constant_coeffs) = self.0.clone().split_off_first();
        let new_constant = constant_coeff * constant;
        let new_non_constant = non_constant_coeffs
            .iter()
            .map(|coeff| coeff * constant)
            .collect::<Vec<_>>();
        let mut new_coeffs = Vec1::new(new_constant);
        new_coeffs.extend(new_non_constant);
        UnivariatePolynomial(new_coeffs)
    }

    pub fn monic(&self) -> Self {
        let last = self.0.last();
        if last.is_zero() {
            self.clone() // 0多項式はそのまま返す
        } else {
            self.mul_constant(&(BigRational::one() / last))
        }
    }

    pub fn substitute(&self, value: &BigRational) -> BigRational {
        self.0
            .iter()
            .rev()
            .fold(BigRational::zero(), |acc, coeff| acc * value + coeff)
    }

    pub fn leading_coeff(&self) -> BigRational {
        self.0.last().clone()
    }

    pub fn get_coeffs(&self) -> &Vec1<BigRational> {
        &self.0
    }
}

impl std::ops::Add for UnivariatePolynomial {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let (constant, non_constant) = self.0.split_off_first();
        let (rhs_constant, rhs_non_constant) = rhs.0.split_off_first();
        let mut sum = Vec1::new(constant + rhs_constant);
        let max_len = usize::max(non_constant.len(), rhs_non_constant.len());
        for i in 0..max_len {
            let coeff1 = non_constant
                .get(i)
                .cloned()
                .unwrap_or_else(BigRational::zero);
            let coeff2 = rhs_non_constant
                .get(i)
                .cloned()
                .unwrap_or_else(BigRational::zero);
            sum.push(coeff1 + coeff2);
        }
        let mut ans = UnivariatePolynomial(sum);
        ans.clean();
        ans
    }
}

impl std::ops::Neg for UnivariatePolynomial {
    type Output = Self;

    fn neg(self) -> Self {
        let (constant, non_constant) = self.0.split_off_first();
        let mut neg_coeffs = Vec1::new(-constant);
        for coeff in non_constant {
            neg_coeffs.push(-coeff);
        }
        UnivariatePolynomial(neg_coeffs)
    }
}

impl std::ops::Sub for UnivariatePolynomial {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

/// 多項式同士の割り算を行って余りを求める
/// ただしcleanがされていることが前提
fn remainder(f: &Vec1<BigRational>, g: &Vec1<BigRational>) -> Option<Vec1<BigRational>> {
    let f_degree = f.len_nonzero().get() - 1;
    let g_degree = g.len_nonzero().get() - 1;
    let lt_f = f.last();
    let lt_g = g.last();
    if g_degree == 0 {
        if lt_g.is_zero() {
            None // 0で割ることはできないためNoneを返す
        } else {
            Some(Vec1::new(BigRational::zero())) // 定数で割るとあまりは常に0になる
        }
    } else if f_degree < g_degree {
        Some(f.clone()) // fの方が次数が小さい場合、割り切れないためfが余りになる
    } else {
        let mut ans = f.clone();
        // まずansの最高次の項を消す
        let leading_coeff = lt_f / lt_g;
        ans[f_degree] = BigRational::zero();
        for i in 0..g_degree {
            let g_coeff = g[i].clone();
            ans[f_degree - g_degree + i] -= leading_coeff.clone() * g_coeff;
        }
        // 末尾の0を削除する
        clean(&mut ans);
        remainder(&ans, g)
    }
}

pub fn uni_poly_remainder(
    f: &UnivariatePolynomial,
    g: &UnivariatePolynomial,
) -> Option<UnivariatePolynomial> {
    remainder(&f.0, &g.0).map(UnivariatePolynomial)
}

fn div(f: &Vec1<BigRational>, g: &Vec1<BigRational>) -> Option<Vec1<BigRational>> {
    let f_degree = f.len_nonzero().get() - 1;
    let g_degree = g.len_nonzero().get() - 1;
    let lt_f = f.last();
    let lt_g = g.last();
    if lt_g.is_zero() {
        None // 0で割ることはできないためNoneを返す
    } else {
        match f_degree.cmp(&g_degree) {
            std::cmp::Ordering::Less => Some(Vec1::new(BigRational::zero())), // fの方が次数が小さい場合、商は常に0になる
            std::cmp::Ordering::Equal => Some(Vec1::new(lt_f / lt_g)), // 同じ次数の場合、商は定数になる
            std::cmp::Ordering::Greater => {
                let res = lt_f / lt_g;
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

pub fn uni_poly_div(
    f: &UnivariatePolynomial,
    g: &UnivariatePolynomial,
) -> Option<UnivariatePolynomial> {
    div(&f.0, &g.0).map(UnivariatePolynomial)
}

pub fn uni_poly_derivative(poly: &UnivariatePolynomial) -> UnivariatePolynomial {
    let mut derivative_coeffs = vec![];
    for (i, coeff) in poly.0.iter().enumerate().skip(1) {
        derivative_coeffs.push(coeff * BigRational::from_integer((i as i64).into()));
    }
    if let Ok(derivative_coeffs) = Vec1::try_from_vec(derivative_coeffs) {
        UnivariatePolynomial(derivative_coeffs)
    } else {
        UnivariatePolynomial(Vec1::new(BigRational::zero()))
    }
}

/// PolynomialからUnivariatePolynomialを作る関数
pub fn polynomial_to_univariate(poly: &Polynomial) -> Result<UnivariatePolynomial> {
    let coeffs =
    poly.raw_iter().map(|(exp, coeff)| {
        match exp.as_slice() {
            [] => Ok((0, coeff.clone())),
            [ind] => Ok((*ind as usize, coeff.clone())),
            _ => Err(color_eyre::eyre::eyre!(
                "Polynomial contains terms with more than one variable, which cannot be converted to UnivariatePolynomial"
            )),
        }
    }).collect::<Result<Vec<_>>>()?;

    let max_ind = coeffs.iter().map(|(ind, _)| *ind).max().unwrap_or(0);
    let mut univariate_coeffs = Vec1::new(BigRational::zero());
    univariate_coeffs.extend(vec![BigRational::zero(); max_ind]);
    for (ind, coeff) in coeffs {
        univariate_coeffs[ind] = coeff;
    }
    Ok(UnivariatePolynomial(univariate_coeffs))
}

impl std::fmt::Display for UnivariatePolynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let degree = self.degree();
        for (i, coeff) in self.0.iter().enumerate().rev() {
            let exp = if i == 0 {
                String::new()
            } else if i == 1 {
                "x".to_string()
            } else {
                format!("x^{}", i)
            };
            let (coeff, is_positive) = if coeff.is_zero() {
                continue;
            } else if coeff == &BigRational::one() && i != 0 {
                (String::new(), true)
            } else if coeff > &BigRational::zero() {
                (format!("{}", coeff), true)
            } else {
                (format!("{}", -coeff), false)
            };
            if i == degree {
                write!(f, "{}{}", coeff, exp)?;
            } else if is_positive {
                write!(f, " + {}{} ", coeff, exp)?;
            } else {
                write!(f, " - {}{} ", coeff, exp)?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use vec1::vec1;

    #[test]
    fn test_remainder() {
        let f = UnivariatePolynomial::new(vec1![
            BigRational::from_integer(9.into()),
            BigRational::from_integer(24.into()),
            BigRational::from_integer(31.into()),
            BigRational::from_integer(23.into()),
            BigRational::from_integer(8.into()),
            BigRational::from_integer(1.into())
        ]);
        let g = UnivariatePolynomial::new(vec1![
            BigRational::from_integer(4.into()),
            BigRational::from_integer(8.into()),
            BigRational::from_integer(5.into()),
            BigRational::from_integer(1.into())
        ]);
        let r = uni_poly_remainder(&f, &g).unwrap();
        let expected = UnivariatePolynomial::new(vec1![
            BigRational::from_integer(9.into()),
            BigRational::from_integer(12.into()),
            BigRational::from_integer(3.into())
        ]);
        assert_eq!(r, expected);
        let r2 = uni_poly_remainder(&g, &r).unwrap();
        assert_eq!(
            r2,
            UnivariatePolynomial::new(vec1![
                BigRational::from_integer(1.into()),
                BigRational::from_integer(1.into())
            ])
        );
        let r3 = uni_poly_remainder(&r, &r2).unwrap();
        assert_eq!(
            r3,
            UnivariatePolynomial::new(vec1![BigRational::from_integer(0.into())])
        );
    }
}
