use color_eyre::Result;
use lalrpop_util::lalrpop_mod;

lalrpop_mod!(
    pub grammar,
    "/parser/grammar.rs"
);

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Expr {
    Num(i64),
    Var(String),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, usize),
    Neg(Box<Expr>),
}

/// Parse an input string into `ast::Expr` using the LALRPOP-generated parser.
pub fn parse_str(input: &str) -> Result<Expr> {
    grammar::ExprParser::new()
        .parse(input)
        .map_err(|e| color_eyre::eyre::eyre!("parse error: {:?}", e))
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::Expr::*;

    #[test]
    fn test_parse_num() {
        let expr = parse_str("-(x^2 + 1)^3 + 42").unwrap();
        let ans = Add(
            Box::new(Neg(Box::new(Pow(
                Box::new(Add(
                    Box::new(Pow(Box::new(Var("x".to_string())), 2)),
                    Box::new(Num(1)),
                )),
                3,
            )))),
            Box::new(Num(42)),
        );
        assert_eq!(expr, ans);

        let expr = parse_str("-x^2 + 1").unwrap();
        let ans = Add(
            Box::new(Neg(Box::new(Pow(Box::new(Var("x".to_string())), 2)))),
            Box::new(Num(1)),
        );
        assert_eq!(expr, ans);
    }
}
