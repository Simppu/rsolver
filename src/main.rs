use std::collections::HashMap;

use rsolver::tokenizer::Tokenizer;
use rsolver::parser::{ Parser};
use rsolver::Expr;




fn main() {
    // Example usage:
    let i = "a+b+a";
    let mut to = Tokenizer::new(i);

    let r = to.tokenize();

    println!("{:?}", r);
    let expr = Expr::Add(vec![
        Expr::Symbol("a".to_string()),
        Expr::Symbol("b".to_string()),
        Expr::Symbol("a".to_string()),
    ]);
    let mut pars = Parser::new(r);

    let e = pars.parse_expression().unwrap();

    println!("{:?}", e.simplify());
    println!("{:?}", expr.simplify());
}