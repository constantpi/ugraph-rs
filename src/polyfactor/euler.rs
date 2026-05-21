use num_traits::{One, Zero};

/// a, bが与えられたときにax+by=gcd(a, b)hを満たすx, yを求める関数
/// 足し算、引き算、掛け算、割り算が定義された環に対して定義できる
pub fn extended_gcd<T>(a: T, b: T) -> (T, T)
where
    T: Clone
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + std::ops::Mul<Output = T>
        + std::ops::Div<Output = T>
        + std::ops::Rem<Output = T>
        + One
        + Zero,
{
    if b.is_zero() {
        (T::one(), T::zero())
    } else {
        let (x1, y1) = extended_gcd(b.clone(), a.clone() % b.clone());
        let x = y1.clone();
        let y = x1 - (a / b) * y1;
        (x, y)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_extended_gcd() {
        let a = 21;
        let b = 30;
        let (x, y) = extended_gcd(a, b);
        println!("{} * {} + {} * {} = {}", a, x, b, y, a * x + b * y);
        assert_eq!(a * x + b * y, 3);
    }
}
