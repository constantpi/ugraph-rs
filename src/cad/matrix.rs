use color_eyre::Result;
use itertools::Itertools;
use num_traits::Zero;
use std::collections::HashSet;
use std::num::NonZeroUsize;
use std::ops::{Add, Mul, Neg};
use vec1::Vec1;

/// Square matrix with elements of type `R`.
#[derive(Clone, Debug)]
pub struct Matrix<R> {
    n: NonZeroUsize,
    data: Vec1<Vec1<R>>,
}

impl<R> Matrix<R>
where
    R: Clone + Zero + Add<Output = R> + Mul<Output = R> + Neg<Output = R> + PartialEq,
{
    /// Create an n x n zero matrix.
    pub fn zero(n: usize) -> Result<Self> {
        let n = NonZeroUsize::new(n)
            .ok_or_else(|| color_eyre::eyre::eyre!("matrix size must be positive"))?;
        let row = Vec1::try_from_vec(vec![R::zero(); n.get()])?;
        Ok(Self {
            n,
            data: Vec1::try_from_vec(vec![row; n.get()]).unwrap(),
        })
    }

    /// Get size
    fn size(&self) -> usize {
        self.n.get()
    }

    /// Indexing (immutable)
    fn get(&self, i: usize, j: usize) -> &R {
        &self.data[i][j]
    }

    pub fn set(&mut self, i: usize, j: usize, value: R) {
        self.data[i][j] = value;
    }

    /// 行列式を計算する
    pub fn determinant(&self, neg_list: &HashSet<Vec<usize>>) -> R {
        let mut det = R::zero();
        for perm in (0..self.size()).permutations(self.size()) {
            let is_neg = neg_list.contains(&perm);
            let perm = Vec1::try_from_vec(perm).unwrap();
            let (first, rest) = perm.split_off_first();
            let mut product = self.get(0, first).clone();
            for (i, j) in rest.into_iter().enumerate() {
                product = product * self.get(i + 1, j).clone();
            }
            if is_neg {
                product = -product;
            }
            det = det + product;
        }
        det
    }

    pub fn get_data(&self) -> Vec<Vec<R>> {
        self.data
            .clone()
            .into_iter()
            .map(|row| row.into_vec())
            .collect()
    }
}

/// permの偶奇を判定する関数
fn is_odd_permutation(perm: &[usize]) -> bool {
    let mut visited = vec![false; perm.len()];
    let mut odd = false;
    for i in 0..perm.len() {
        if !visited[i] {
            let mut cycle_length = 0;
            let mut j = i;
            while !visited[j] {
                visited[j] = true;
                j = perm[j];
                cycle_length += 1;
            }
            if cycle_length % 2 == 0 {
                odd = !odd;
            }
        }
    }
    odd
}

pub fn generate_neg_list(n: usize) -> HashSet<Vec<usize>> {
    (0..n)
        .permutations(n)
        .filter(|perm| is_odd_permutation(perm))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determinant() -> Result<()> {
        let mat_vec = vec![vec![8, 5, -3], vec![0, 7, 3], vec![-2, -1, 1]];
        let mut m = Matrix::zero(3)?;
        for i in 0..3 {
            for j in 0..3 {
                m.set(i, j, mat_vec[i][j]);
            }
        }
        let neg_list = generate_neg_list(3);
        assert_eq!(m.determinant(&neg_list), 8);
        Ok(())
    }
}
