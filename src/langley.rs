use color_eyre::{Result, eyre::Ok};

// ラングレーの問題を多項式制約で表現するためのコード

/// cos^2(theta)を表す多項式を生成する関数
fn cos_squared(theta: i32, var: &str) -> Result<String> {
    // cos^2(theta)=(1+cos(2*theta))/2を利用して、thetaに応じた多項式を返す
    let theta_2 = (theta * 2) % 360;
    let theta_2 = if theta_2 < 0 { theta_2 + 360 } else { theta_2 };
    // 0~180度の範囲だけ考えれば十分なので、theta_2を0~180度の範囲に変換する
    let theta_2 = if theta_2 > 180 {
        360 - theta_2
    } else {
        theta_2
    };
    // (2 * var - 1)をひとかたまりとする。
    let var_squared = format!("(2 * {} - 1)", var);
    // theta_2が20度の倍数である場合にのみ対応する
    // cos3θ=4cos^3θ-3cosθを利用して、theta_2に応じた多項式を返す
    match theta_2 {
        0 => Ok(format!("{var_squared}-1")),
        20 | 100 | 140 => Ok(format!("8*{var_squared}^3-6*{var_squared}-1")),
        40 | 80 | 160 => Ok(format!("8*{var_squared}^3-6*{var_squared}+1")),
        60 => Ok(format!("2*{var_squared}-1")),
        120 => Ok(format!("2*{var_squared}+1")),
        180 => Ok(format!("{var_squared}+1")),
        _ => Err(color_eyre::eyre::eyre!(
            "theta must be a multiple of 20 degrees, but got {theta}"
        )),
    }
}

/// 内積の制約を表す多項式を生成する関数
fn inner_product_constraint(
    a: (&str, &str),
    o: (&str, &str),
    b: (&str, &str),
    cos2_theta: &str,
) -> String {
    // まずはベクトルを定義する
    let a_x = format!("({} - {})", a.0, o.0);
    let a_y = format!("({} - {})", a.1, o.1);
    let b_x = format!("({} - {})", b.0, o.0);
    let b_y = format!("({} - {})", b.1, o.1);
    //  (内積)^2 - cos^2(theta) * (|a|^2 * |b|^2) = 0を表す多項式を生成する
    let inner_product_squared = format!("({a_x}*{b_x}+{a_y}*{b_y})^2");
    let a_squared = format!("({a_x}^2 + {a_y}^2)");
    let b_squared = format!("({b_x}^2 + {b_y}^2)");
    format!("{inner_product_squared} - {cos2_theta} * {a_squared} * {b_squared}")
}

pub fn generate_langley_polynomials() -> Result<Vec<String>> {
    // 左上から反時計回りに、点A, B, C, Dを配置する
    let a = ("a_x", "a_y");
    let b = ("b_x", "b_y");
    let c = ("c_x", "c_y");
    let d = ("d_x", "d_y");

    let cos2_20 = "cos2_20";
    let cos2_30 = "cos2_30";
    let cos2_50 = "cos2_50";
    let cos2_60 = "cos2_60";

    let mut constraints = Vec::new();
    constraints.push(cos_squared(20, cos2_20)?);
    constraints.push(cos_squared(30, cos2_30)?);
    constraints.push(cos_squared(50, cos2_50)?);
    constraints.push(cos_squared(60, cos2_60)?);

    // 内積の制約を追加する
    // 角CBD = 60度
    constraints.push(inner_product_constraint(c, b, d, cos2_60));
    // 角DBA = 20度
    constraints.push(inner_product_constraint(d, b, a, cos2_20));
    // 角DCA = 30度
    constraints.push(inner_product_constraint(d, c, a, cos2_30));
    // 角ACB = 50度
    constraints.push(inner_product_constraint(a, c, b, cos2_50));

    // 最後に基準点を固定する
    // 点Bを原点に固定する
    constraints.push(format!("{}", b.0));
    constraints.push(format!("{}", b.1));
    // 点Cをx軸上に固定する
    constraints.push(format!("{}-1", c.0));
    constraints.push(format!("{}", c.1));
    Ok(constraints)
}
