use std::ops::{Add, AddAssign, Div, DivAssign, Mul, Neg, Sub, SubAssign};

/// 素数上の体
#[derive(Clone, Debug, PartialEq, Eq, Copy)]
pub struct PrimeField {
    n: usize,
    p: usize,
}

impl PrimeField {
    /// 素数pで初期化
    pub fn new(n: usize, p: usize) -> Self {
        if !is_prime(p) {
            panic!("p must be prime");
        }
        Self { n: n % p, p }
    }
    fn clear(&mut self) {
        self.n %= self.p;
    }
    pub fn get_prime(&self) -> usize {
        self.p
    }
    pub fn is_zero(&self) -> bool {
        self.n == 0
    }
    fn pow(&self, exp: usize) -> Self {
        match exp {
            0 => Self::new(1, self.p),
            1 => *self,
            _ => {
                let half = self.pow(exp / 2);
                let result = half * half;
                if exp % 2 == 1 { result * *self } else { result }
            }
        }
    }
    pub fn zero(p: usize) -> Self {
        Self::new(0, p)
    }
    pub fn one(p: usize) -> Self {
        Self::new(1, p)
    }
}

impl Add for PrimeField {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.p != rhs.p {
            panic!("Cannot add elements from different fields");
        }
        let mut result = Self::new((self.n + rhs.n) % self.p, self.p);
        result.clear();
        result
    }
}

impl AddAssign for PrimeField {
    fn add_assign(&mut self, rhs: Self) {
        if self.p != rhs.p {
            panic!("Cannot add elements from different fields");
        }
        self.n = (self.n + rhs.n) % self.p;
        self.clear();
    }
}

impl Neg for PrimeField {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self::new((self.p - self.n) % self.p, self.p)
    }
}

impl Sub for PrimeField {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

impl SubAssign for PrimeField {
    fn sub_assign(&mut self, rhs: Self) {
        *self += -rhs;
    }
}

impl Mul for PrimeField {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        if self.p != rhs.p {
            panic!("Cannot multiply elements from different fields");
        }
        Self::new((self.n * rhs.n) % self.p, self.p)
    }
}

impl Div for PrimeField {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        if self.p != rhs.p {
            panic!("Cannot divide elements from different fields");
        }
        self * rhs.pow(self.p - 2)
    }
}

pub fn is_prime(n: usize) -> bool {
    if n <= 1 {
        false
    } else {
        let mut i = 2;
        loop {
            if i * i > n {
                break true;
            }
            if n % i == 0 {
                break false;
            }
            i += 1;
        }
    }
}

/// 素数を生成するiterator
pub struct PrimeIter {
    current: usize,
}
impl PrimeIter {
    pub fn new() -> Self {
        Self { current: 1 }
    }
}

impl Iterator for PrimeIter {
    type Item = usize;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            self.current += 1;
            if is_prime(self.current) {
                return Some(self.current);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_prime_iter() {
        let mut prime_iter = PrimeIter::new();
        let primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29];
        for &p in primes.iter() {
            assert_eq!(prime_iter.next(), Some(p));
        }
    }
}
