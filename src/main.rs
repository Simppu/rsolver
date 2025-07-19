use std::collections::HashMap;

use rsolver::tokenizer::Tokenizer;
use rsolver::parser::{ Parser};
use rsolver::{evaluate_expr, evaluate_expr1, Expr};




fn main() {
    // Example usage:
    let i = "1/((1-1)) + a";
    let mut to = Tokenizer::new(i);

    let r = to.tokenize();

    
    
    let mut pars = Parser::new(r);

    let e = pars.parse_expression().unwrap();
    
    println!("{}", evaluate_expr1(e.simplify()));
}