use color_eyre::Result;

use super::{RelOp, ast_to_polynomial, parse_str};
use crate::polynomial::Polynomial;

/// ファイルを読み込んで、そこに書いてある式をPolynomialに変換する
pub fn read_file_to_polynomial(path: &str) -> Result<Vec<(Polynomial, RelOp)>> {
    let content = std::fs::read_to_string(path)?;
    let mut expr_list = Vec::new();
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let expr = parse_str(line)?;
        expr_list.push(expr);
    }
    ast_to_polynomial(&expr_list)
}
