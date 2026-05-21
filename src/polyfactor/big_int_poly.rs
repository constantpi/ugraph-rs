use num::{BigInt, FromPrimitive, One, Zero};
use vec1::{Vec1, vec1};

fn clean(coefficients: &mut Vec1<BigInt>) {
    // popが成功すれば残りもVec1であるため続行する
    if coefficients.last().is_zero() && coefficients.pop().is_ok() {
        clean(coefficients);
    }
}

/// 一変数多項式を表す構造体
/// i番目はx^iの係数を表す
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BigIntPoly(Vec1<BigInt>);

impl BigIntPoly {
    pub fn new(coefficients: Vec1<BigInt>) -> Self {
        let mut ans = BigIntPoly(coefficients);
        ans.clean();
        ans
    }

    pub fn degree(&self) -> usize {
        self.0.len_nonzero().get() - 1
    }

    pub fn iter(&self) -> impl Iterator<Item = &BigInt> {
        self.0.iter()
    }

    /// 末尾の0を削除する関数
    pub fn clean(&mut self) {
        clean(&mut self.0);
    }

    pub fn is_zero(&self) -> bool {
        self.0.len_nonzero().get() == 1 && self.0.first().is_zero()
    }

    pub fn is_one(&self) -> bool {
        self.0.len_nonzero().get() == 1 && self.0.first().is_one()
    }

    pub fn mul_constant(&self, constant: &BigInt) -> Self {
        let (constant_coeff, non_constant_coeffs) = self.0.clone().split_off_first();
        let new_constant = constant_coeff * constant;
        let new_non_constant = non_constant_coeffs
            .iter()
            .map(|coeff| coeff * constant)
            .collect::<Vec<_>>();
        let mut new_coeffs = Vec1::new(new_constant);
        new_coeffs.extend(new_non_constant);
        BigIntPoly(new_coeffs)
    }

    pub fn mod_integer(&self, modulus: &BigInt) -> Self {
        let (rest, last) = self.0.clone().split_off_last();
        let new_rest = rest.iter().map(|coeff| coeff % modulus).collect::<Vec<_>>();
        let new_last = last % modulus;
        let new_coeffs = Vec1::from_vec_push(new_rest, new_last);
        BigIntPoly(new_coeffs)
    }
}

impl std::ops::Add for BigIntPoly {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        let (constant, non_constant) = self.0.split_off_first();
        let (rhs_constant, rhs_non_constant) = rhs.0.split_off_first();
        let mut sum = Vec1::new(constant + rhs_constant);
        let max_len = usize::max(non_constant.len(), rhs_non_constant.len());
        for i in 0..max_len {
            let coeff1 = non_constant.get(i).cloned().unwrap_or_else(BigInt::zero);
            let coeff2 = rhs_non_constant
                .get(i)
                .cloned()
                .unwrap_or_else(BigInt::zero);
            sum.push(coeff1 + coeff2);
        }
        let mut ans = BigIntPoly(sum);
        ans.clean();
        ans
    }
}

impl std::ops::Neg for BigIntPoly {
    type Output = Self;

    fn neg(self) -> Self {
        let (constant, non_constant) = self.0.split_off_first();
        let mut neg_coeffs = Vec1::new(-constant);
        for coeff in non_constant {
            neg_coeffs.push(-coeff);
        }
        BigIntPoly(neg_coeffs)
    }
}

impl std::ops::Sub for BigIntPoly {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        self + (-rhs)
    }
}

impl std::ops::Mul for BigIntPoly {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        // let mut result_coeffs = vec1![BigInt::zero(); self.degree() + rhs.degree() + 1];
        let mut result_coeffs = Vec1::new(BigInt::zero());
        result_coeffs.extend(vec![BigInt::zero(); self.degree() + rhs.degree()]);
        for (i, coeff1) in self.0.iter().enumerate() {
            for (j, coeff2) in rhs.0.iter().enumerate() {
                result_coeffs[i + j] += coeff1 * coeff2;
            }
        }
        BigIntPoly(result_coeffs)
    }
}

impl std::ops::Div for BigIntPoly {
    type Output = Self;

    fn div(self, rhs: Self) -> Self {
        uni_poly_div(&self, &rhs).expect("割り算に失敗しました")
    }
}

impl std::ops::Rem for BigIntPoly {
    type Output = Self;

    fn rem(self, rhs: Self) -> Self {
        uni_poly_remainder(&self, &rhs).expect("割り算に失敗しました")
    }
}

impl num_traits::Zero for BigIntPoly {
    fn zero() -> Self {
        BigIntPoly(Vec1::new(BigInt::zero()))
    }

    fn is_zero(&self) -> bool {
        self.is_zero()
    }
}

impl num_traits::One for BigIntPoly {
    fn one() -> Self {
        BigIntPoly(Vec1::new(BigInt::one()))
    }

    fn is_one(&self) -> bool {
        self.is_one()
    }
}

/// 多項式同士の割り算を行って余りを求める
/// ただしcleanがされていることが前提
/// またgの最高次の係数が1であることが前提
fn remainder(f: &Vec1<BigInt>, g: &Vec1<BigInt>) -> Option<Vec1<BigInt>> {
    let f_degree = f.len_nonzero().get() - 1;
    let g_degree = g.len_nonzero().get() - 1;
    let lt_f = f.last();
    let lt_g = g.last();
    if !lt_g.is_one() {
        None // gの最高次の係数が1でない場合は割り算できないためNoneを返す
    } else {
        if g_degree == 0 {
            Some(Vec1::new(BigInt::zero())) // 1で割るとあまりは常に0になる
        } else if f_degree < g_degree {
            Some(f.clone()) // fの方が次数が小さい場合、割り切れないためfが余りになる
        } else {
            let mut ans = f.clone();
            // まずansの最高次の項を消す
            let leading_coeff = lt_f;
            ans[f_degree] = BigInt::zero();
            for i in 0..g_degree {
                let g_coeff = g[i].clone();
                ans[f_degree - g_degree + i] -= leading_coeff.clone() * g_coeff;
            }
            // 末尾の0を削除する
            clean(&mut ans);
            remainder(&ans, g)
        }
    }
}

pub fn uni_poly_remainder(f: &BigIntPoly, g: &BigIntPoly) -> Option<BigIntPoly> {
    remainder(&f.0, &g.0).map(BigIntPoly)
}

pub fn uni_poly_division(f: &BigIntPoly, g: &BigIntPoly) -> Option<BigIntPoly> {
    div(&f.0, &g.0).map(BigIntPoly)
}

fn div(f: &Vec1<BigInt>, g: &Vec1<BigInt>) -> Option<Vec1<BigInt>> {
    let f_degree = f.len_nonzero().get() - 1;
    let g_degree = g.len_nonzero().get() - 1;
    let lt_f = f.last();
    let lt_g = g.last();
    if !lt_g.is_one() {
        // gの最高次の係数が1でない場合は割り算できないため
        None
    } else {
        match f_degree.cmp(&g_degree) {
            std::cmp::Ordering::Less => Some(Vec1::new(BigInt::zero())), // fの方が次数が小さい場合、商は常に0になる
            std::cmp::Ordering::Equal => Some(Vec1::new(lt_f.clone())), // 同じ次数の場合、商は定数になる
            std::cmp::Ordering::Greater => {
                let res = lt_f;
                let mut f = f.clone();
                f.pop().unwrap(); // 最高次の項を削除。fの次数はg以上なので必ず成功する
                for i in 0..g_degree {
                    let g_coeff = g[i].clone();
                    f[f_degree - g_degree + i] -= res.clone() * g_coeff;
                }
                let mut ans = div(&f, g)?;
                ans.push(res.clone());
                Some(ans)
            }
        }
    }
}

pub fn uni_poly_div(f: &BigIntPoly, g: &BigIntPoly) -> Option<BigIntPoly> {
    div(&f.0, &g.0).map(BigIntPoly)
}

pub fn uni_poly_derivative(poly: &BigIntPoly) -> BigIntPoly {
    let mut derivative_coeffs = vec![];
    for (i, coeff) in poly.0.iter().enumerate().skip(1) {
        derivative_coeffs.push(coeff * BigInt::from_usize(i).unwrap());
    }
    if let Some(derivative_coeffs) = Vec1::try_from_vec(derivative_coeffs).ok() {
        BigIntPoly(derivative_coeffs)
    } else {
        BigIntPoly(Vec1::new(BigInt::zero()))
    }
}
