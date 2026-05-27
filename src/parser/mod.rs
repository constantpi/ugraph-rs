mod ast;
mod gen_poly;
mod read_file;

use ast::{Expr, parse_str};
use gen_poly::ast_to_polynomial;
pub use read_file::read_file_to_polynomial;
