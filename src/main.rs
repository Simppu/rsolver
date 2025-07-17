use std::collections::HashMap;

use rsolver::tokenizer::Tokenizer;
use rsolver::parser::{ Parser};




fn main() {
    
    let i = "a+b+a";
    let mut to = Tokenizer::new(i);

    let r = to.tokenize();

    println!("{:?}", r);

    let mut pars = Parser::new(r);

    let e = pars.parse_expression().unwrap();

    println!("{:?}", e);
    println!("{:?}", e.simplify());
}