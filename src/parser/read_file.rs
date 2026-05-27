use color_eyre::Result;

use super::{ast_to_polynomial, parse_str};
use crate::polynomial::Polynomial;

/// ファイルを読み込んで、そこに書いてある式をPolynomialに変換する
pub fn read_file_to_polynomial(path: &str) -> Result<Vec<Polynomial>> {
    let content = std::fs::read_to_string(path)?;
    let mut polynomials = Vec::new();
    for line in content.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let expr = parse_str(line)?;
        let poly = ast_to_polynomial(&expr)?;
        polynomials.push(poly);
    }
    Ok(polynomials)
}
