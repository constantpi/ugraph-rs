use color_eyre::Result;
use lalrpop_util::lalrpop_mod;

lalrpop_mod!(
    pub grammar,
    "/parser/grammar.rs"
);

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum RelOp {
    Eq,
    Neq,
    Lt,
    Gt,
    Leq,
    Geq,
}

impl std::fmt::Display for RelOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            RelOp::Eq => "=",
            RelOp::Neq => "!=",
            RelOp::Lt => "<",
            RelOp::Gt => ">",
            RelOp::Leq => "<=",
            RelOp::Geq => ">=",
        };
        write!(f, "{}", s)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Inequality {
    left: Expr,
    op: RelOp,
    right: Expr,
}

impl Inequality {
    pub fn get_left(&self) -> &Expr {
        &self.left
    }
    pub fn get_op(&self) -> RelOp {
        self.op
    }
    pub fn get_right(&self) -> &Expr {
        &self.right
    }
}

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
pub fn parse_str(input: &str) -> Result<Inequality> {
    grammar::InequalityParser::new()
        .parse(input)
        .map_err(|e| color_eyre::eyre::eyre!("parse error: {:?}", e))
}

#[cfg(test)]
mod tests {
    use super::*;

    use super::Expr::*;

    #[test]
    fn test_parse_num() {
        let expr = parse_str("-(x^2 + 1)^3 + 42 = 0").unwrap();
        let left = Add(
            Box::new(Neg(Box::new(Pow(
                Box::new(Add(
                    Box::new(Pow(Box::new(Var("x".to_string())), 2)),
                    Box::new(Num(1)),
                )),
                3,
            )))),
            Box::new(Num(42)),
        );
        let right = Num(0);
        let ans = Inequality {
            left,
            op: RelOp::Eq,
            right,
        };
        assert_eq!(expr, ans);

        let expr = parse_str("-x^2 + 1 = 0").unwrap();
        let left = Add(
            Box::new(Neg(Box::new(Pow(Box::new(Var("x".to_string())), 2)))),
            Box::new(Num(1)),
        );
        let right = Num(0);
        let ans = Inequality {
            left,
            op: RelOp::Eq,
            right,
        };
        assert_eq!(expr, ans);
    }
}
