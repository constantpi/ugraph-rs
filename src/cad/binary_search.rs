use itertools::Itertools;
use num::{BigRational, Signed};

use super::{Range, Root, UnivariatePolynomial};

/// 区間ニュートン法を用いて、根を含む区間の精度を高める関数。
fn interval_newton(
    f: &UnivariatePolynomial,
    lower: &BigRational,
    upper: &BigRational,
) -> Option<(BigRational, BigRational)> {
    // ニュートン法の更新式は x_midを有理数として、x_mid - f(x_mid) / f'(Interval)を計算することになる。
    let (is_positive, abs_lower, abs_upper) = if lower.is_positive() {
        (true, lower.clone(), upper.clone())
    } else if upper.is_negative() {
        (false, -upper.clone(), -lower.clone())
    } else {
        // 0を含む場合はニュートン法は適用できない
        return None;
    };
    // 微分値の区間評価を計算する
    let mut pow_positive = true;
    let mut der_lower = BigRational::from_integer(0.into());
    let mut der_upper = BigRational::from_integer(0.into());
    let mut pow_lower = BigRational::from_integer(1.into());
    let mut pow_upper = BigRational::from_integer(1.into());
    for (i, coeff) in f.iter().enumerate().skip(1) {
        // f'(x) = sum(i * coeff[i] * x^(i-1))なので、i=1のときは定数項になる。i=0のときは微分して0になるので無視していい。
        let coeff_abs = coeff.abs();
        if pow_positive == (coeff.is_positive()) {
            der_lower +=
                BigRational::from_integer(i.into()) * coeff_abs.clone() * pow_lower.clone();
            der_upper +=
                BigRational::from_integer(i.into()) * coeff_abs.clone() * pow_upper.clone();
        } else {
            der_lower -=
                BigRational::from_integer(i.into()) * coeff_abs.clone() * pow_upper.clone();
            der_upper -=
                BigRational::from_integer(i.into()) * coeff_abs.clone() * pow_lower.clone();
        }
        pow_lower *= abs_lower.clone();
        pow_upper *= abs_upper.clone();
        if !is_positive {
            pow_positive = !pow_positive;
        }
    }
    let (der_is_positive, der_abs_lower, der_abs_upper) = if der_lower.is_positive() {
        (true, der_lower.clone(), der_upper.clone())
    } else if der_upper.is_negative() {
        (false, -der_upper.clone(), -der_lower.clone())
    } else {
        // 微分値が0を含む場合はニュートン法は適用できない
        return None;
    };

    let midpoint = (lower.clone() + upper.clone()) / BigRational::from_integer(2.into());
    let eval_midpoint = f.substitute(&midpoint);
    let eval_midpoint_abs = eval_midpoint.abs();
    let delta_abs_lower = eval_midpoint_abs.clone() / der_abs_upper.clone();
    let delta_abs_upper = eval_midpoint_abs.clone() / der_abs_lower.clone();
    let delta_is_positive = eval_midpoint.is_positive() == der_is_positive;
    let (new_lower, new_upper) = if delta_is_positive {
        (
            midpoint.clone() - delta_abs_upper.clone(),
            midpoint.clone() - delta_abs_lower.clone(),
        )
    } else {
        (
            midpoint.clone() + delta_abs_lower.clone(),
            midpoint.clone() + delta_abs_upper.clone(),
        )
    };
    Some((new_lower, new_upper))
}

/// 解を含む範囲の精度を高める関数。
/// fは最小多項式であって2次以上であること。
/// lower < upperであり、この間にfの根が1つだけ存在することが保証されていること。
pub fn refine_range(
    f: &UnivariatePolynomial,
    lower: &BigRational,
    upper: &BigRational,
    use_newton: bool,
) -> (BigRational, BigRational) {
    // 順にlower, midpoint, upperを評価して符号が変わるところを探す
    // 評価結果が0になることはありえない。もしそうなら次数が1となるから。
    let midpoint = (lower.clone() + upper.clone()) / BigRational::from_integer(2.into());
    let eval_lower = f.substitute(lower);
    let eval_midpoint = f.substitute(&midpoint);
    let (new_lower, new_upper) = if eval_lower.signum() != eval_midpoint.signum() {
        (lower.clone(), midpoint)
    } else {
        (midpoint, upper.clone())
    };
    if use_newton
        && let Some((newton_lower, newton_upper)) = interval_newton(f, &new_lower, &new_upper)
    {
        // 2つの区間の交わりを取る
        let final_lower = new_lower.clone().max(newton_lower.clone());
        let final_upper = new_upper.clone().min(newton_upper.clone());
        (final_lower, final_upper)
    } else {
        (new_lower, new_upper)
    }
}

/// RootがIntervalのときにその精度を高める
pub fn refine_root(root: &Root, use_newton: bool) -> Option<Root> {
    if let Some((lower, upper)) = root.get_interval() {
        // この時点でlower < upperは保証されている
        // また最小多項式の次数は2以上であることも保証されている
        let poly = root.get_poly();
        let (new_low, new_high) = refine_range(poly, &lower, &upper, use_newton);
        Some(Root::new(poly.clone(), new_low, new_high))
    } else {
        None
    }
}

/// 2つのRootが分離しているかどうかの判定。
fn is_separated(root1: &Root, root2: &Root) -> bool {
    match (root1.get_range(), root2.get_range()) {
        (Range::Exact(r1), Range::Exact(r2)) => r1 != r2, // 根が等しい場合は分離していない
        (Range::Exact(r1), Range::Interval(l2, u2))
        | (Range::Interval(l2, u2), Range::Exact(r1)) => {
            r1 < l2 || r1 > u2 // 根が区間の外にある場合は分離している
        }
        (Range::Interval(l1, u1), Range::Interval(l2, u2)) => {
            u1 < l2 || u2 < l1 // 区間同士が重ならない場合は分離している
        }
    }
}

/// Root集合が与えられる。同じものは含まれないことが保証されている。
/// ここからサンプル点列を作る。サンプル点列はRoot集合のすべてのRootを挟むように選ぶ。
pub fn calc_sample_points(roots: &[Root]) -> Vec<Root> {
    // まずはすべてのRootを完全に被りがないようにする
    let mut distinct_roots = Vec::new();
    for root in roots {
        let mut root = root.clone();
        for other in distinct_roots.iter_mut() {
            // rootとotherが分離していない場合、rootを精度を高めて分離させる
            while !is_separated(&root, other) {
                if let Some(refined) = refine_root(&root, false) {
                    root = refined;
                }
                if let Some(refined) = refine_root(other, false) {
                    *other = refined;
                }
            }
        }
        distinct_roots.push(root);
    }
    // distinct_rootsはすべてのRootが分離していることが保証されている
    // ここから大きさ順にソートする。Exactならその値、Intervalなら下端を基準にする。
    distinct_roots.sort_by_key(|r| r.get_lower_bound());
    // r1, r2があったら、s0をr1の下端-1、s1を(r1の上端+r2の下端)/2、s2をr2の上端+1とする。これでr1とr2を挟む3点ができる。
    let mut sample_points = Vec::new();
    if let Some(first) = distinct_roots.first() {
        let first_lower = first.get_lower_bound();
        sample_points.push(Root::new_rational(
            first_lower - BigRational::from_integer(1.into()),
        ));
        sample_points.push(first.clone());
    } else {
        // 根がない場合は適当に0をサンプル点として追加しておく
        sample_points.push(Root::new_rational(BigRational::from_integer(0.into())));
    }
    for (s, t) in distinct_roots.iter().tuple_windows() {
        let s_upper = s.get_upper_bound();
        let t_lower = t.get_lower_bound();
        let midpoint = (s_upper.clone() + t_lower.clone()) / BigRational::from_integer(2.into());
        sample_points.push(Root::new_rational(midpoint));
        sample_points.push(t.clone());
    }
    if let Some(last) = distinct_roots.last() {
        let last_upper = last.get_upper_bound();
        sample_points.push(Root::new_rational(
            last_upper + BigRational::from_integer(1.into()),
        ));
    }

    sample_points
}
