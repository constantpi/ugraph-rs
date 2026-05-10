use num::{BigRational, Zero};
use vec1::Vec1;

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
