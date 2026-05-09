use num::{BigRational, Zero};
use std::collections::HashMap;

/// n変数の多項式の項の指数を表す構造体
/// 例えば、x^2 * y^3 の場合、Exponentは [2, 3] となる
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct Exponent(Vec<u32>);

/// n変数の多項式を表す構造体
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Polynomial {
    // 項の指数と係数のマップ
    terms: HashMap<Exponent, BigRational>,
}

impl Exponent {
    pub fn new(exponents: Vec<u32>) -> Self {
        Exponent(exponents)
    }

    pub fn sub(&self, other: &Self) -> Option<Self> {
        if self.0.len() != other.0.len() {
            panic!("異なる変数の数の指数は減算できません");
        }
        if self.0.iter().zip(other.0.iter()).any(|(a, b)| a < b) {
            None // 減算結果が負の指数になる場合はNoneを返す
        } else {
            let new_exponents = self
                .0
                .iter()
                .zip(other.0.iter())
                .map(|(a, b)| a - b)
                .collect();
            Some(Exponent(new_exponents))
        }
    }

    pub fn saturated_sub(&self, other: &Self) -> Self {
        if self.0.len() != other.0.len() {
            panic!("異なる変数の数の指数は減算できません");
        }
        let new_exponents = self
            .0
            .iter()
            .zip(other.0.iter())
            .map(|(a, b)| a.saturating_sub(*b))
            .collect();
        Exponent(new_exponents)
    }

    pub fn lcm(&self, other: &Self) -> Self {
        if self.0.len() != other.0.len() {
            panic!("異なる変数の数の指数はLCMを計算できません");
        }
        let new_exponents = self
            .0
            .iter()
            .zip(other.0.iter())
            .map(|(a, b)| std::cmp::max(*a, *b))
            .collect();
        Exponent(new_exponents)
    }

    pub fn sum_degree(&self) -> u32 {
        self.0.iter().sum()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    /// 一番後ろの変数の次数とその変数を除いた指数を返す
    pub fn split_last(&self) -> Option<(u32, Self)> {
        let mut clone = self.clone();
        let last = clone.0.pop()?;
        Some((last, Exponent(clone.0)))
    }

    /// 一番後ろに変数の次数を追加する
    pub fn push(&mut self, exp: u32) {
        self.0.push(exp);
    }
}
// ExponentにAdd トレイトを実装
impl std::ops::Add for Exponent {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        if self.0.len() != other.0.len() {
            panic!("異なる変数の数の指数は加算できません");
        }
        let new_exponents = self
            .0
            .iter()
            .zip(other.0.iter())
            .map(|(a, b)| a + b)
            .collect();
        Exponent(new_exponents)
    }
}

// ExponentにLexicographic Orderを実装（全順序を定義）
impl PartialOrd for Exponent {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Exponent {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let len_cmp = self.0.len().cmp(&other.0.len());
        if len_cmp != std::cmp::Ordering::Equal {
            return len_cmp;
        }
        for (a, b) in self.0.iter().zip(other.0.iter()) {
            match a.cmp(b) {
                std::cmp::Ordering::Less => return std::cmp::Ordering::Less,
                std::cmp::Ordering::Greater => return std::cmp::Ordering::Greater,
                std::cmp::Ordering::Equal => continue,
            }
        }
        std::cmp::Ordering::Equal
    }
}

impl Polynomial {
    pub fn zero() -> Self {
        Polynomial {
            terms: HashMap::new(),
        }
    }

    fn clean(&mut self) {
        self.terms.retain(|_, coeff| !coeff.is_zero());
    }

    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// 多項式に項を追加するメソッド
    pub fn add_term(&mut self, exponent: Exponent, coefficient: BigRational) {
        *self
            .terms
            .entry(exponent)
            .or_insert(BigRational::from_integer(0.into())) += coefficient;
        self.clean();
    }

    /// iteratorを返すメソッド
    /// ただし、項の順序はLexicographic Orderでソートされている必要がある
    pub fn lex_iter(&self) -> impl Iterator<Item = (&Exponent, &BigRational)> {
        let mut terms = self.raw_iter().collect::<Vec<_>>();
        terms.sort_by(|(exp_a, _), (exp_b, _)| exp_b.cmp(exp_a)); // 降順でソート
        terms.into_iter()
    }

    pub fn raw_iter(&self) -> impl Iterator<Item = (&Exponent, &BigRational)> {
        self.terms.iter()
    }

    pub fn mul_term(&self, exponent: Exponent, coefficient: BigRational) -> Self {
        if coefficient.is_zero() {
            Polynomial::zero()
        } else {
            let new_terms = self
                .raw_iter()
                .map(|(exp, coeff)| (exp.clone() + exponent.clone(), coeff * &coefficient))
                .collect::<HashMap<_, _>>();
            Polynomial { terms: new_terms }
        }
    }
    pub fn mul_rational(&self, coefficient: BigRational) -> Self {
        if coefficient.is_zero() {
            Polynomial::zero()
        } else {
            let new_terms = self
                .raw_iter()
                .map(|(exp, coeff)| (exp.clone(), coeff * &coefficient))
                .collect::<HashMap<_, _>>();
            Polynomial { terms: new_terms }
        }
    }

    pub fn get_lt(&self) -> Option<(&Exponent, &BigRational)> {
        self.lex_iter().next()
    }
}

// PolynomialにAdd トレイトを実装
impl std::ops::Add for Polynomial {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        let mut result = self.clone();
        for (exp, coeff) in other.raw_iter() {
            result.add_term(exp.clone(), coeff.clone());
        }
        result
    }
}

// PolynomialにNeg トレイトを実装
impl std::ops::Neg for Polynomial {
    type Output = Self;
    fn neg(self) -> Self {
        self.mul_rational(BigRational::from_integer((-1).into()))
    }
}

// PolynomialにSub トレイトを実装
impl std::ops::Sub for Polynomial {
    type Output = Self;
    fn sub(self, other: Self) -> Self {
        self + (-other)
    }
}

// PolynomialにMul トレイトを実装
impl std::ops::Mul for Polynomial {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        let mut terms = HashMap::new();
        for (exp_a, coeff_a) in self.raw_iter() {
            for (exp_b, coeff_b) in other.raw_iter() {
                let new_exp = exp_a.clone() + exp_b.clone();
                let new_coeff = coeff_a * coeff_b;
                *terms
                    .entry(new_exp)
                    .or_insert(BigRational::from_integer(0.into())) += new_coeff;
            }
        }
        let mut ans = Polynomial { terms };
        ans.clean();
        ans
    }
}

impl std::ops::AddAssign for Polynomial {
    fn add_assign(&mut self, other: Self) {
        for (exp, coeff) in other.raw_iter() {
            self.add_term(exp.clone(), coeff.clone());
        }
    }
}

impl std::fmt::Display for Exponent {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for (variables, exponent) in self.0.iter().enumerate() {
            if *exponent == 1 {
                write!(f, "x{}", variables + 1)?;
            } else if *exponent > 1 {
                write!(f, "x{}^{}", variables + 1, exponent)?;
            }
        }
        Ok(())
    }
}

// Display トレイトを実装
impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let Some(first) = self.lex_iter().next() else {
            return write!(f, "0");
        };
        write!(f, "{}{}", first.1, first.0)?;
        for (exp, coeff) in self.lex_iter().skip(1) {
            if coeff.is_zero() {
                continue;
            }
            if coeff > &BigRational::from_integer(0.into()) {
                write!(f, " + {}{}", coeff, exp)?;
            } else {
                write!(f, " - {}{}", -coeff, exp)?;
            }
        }
        Ok(())
    }
}
