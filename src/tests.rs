#[cfg(test)]
use crate::Expr;
use crate::{parser::Parser, tokenizer::{Rational, Token, Tokenizer}};



#[test]
pub fn rational_ident() {
    let a = Rational::new(5, 1);
    let b = Rational::new(5, 1);

    assert!(a == b);

   

}

#[test]
pub fn rational_add() {
    let a = Rational::new(5, 1);
    let b = Rational::new(12, 7);

    assert_eq!(a + b, Rational::new(47, 7))
}
#[test]
pub fn rational_sub() {
    let a = Rational::new(5, 1);
    let b = Rational::new(12, 7);

    assert_eq!(a - b, Rational::new(23, 7))
}



#[test]
pub fn tokenizer() {
    let input = "5";
    let mut tokenizer = Tokenizer::new(input);
    let a = tokenizer.tokenize();
    let nums = vec![
        Token::Number(Rational::new(5, 1))
    ];
    println!("{:?}", a);
    println!("{:?}", nums);
    println!("{}", nums == a);



}

#[test]
pub fn algebraic_addition() {
    let input = "a + a";
    assert_eq!(format!("{}", expression(input)), "2 * a")
}

#[test]
pub fn algebraic_multiplication() {
    let input = "a * a";
    assert_eq!(format!("{}", expression(input)), "a^2")
}

#[test]
pub fn algebraic_division() {
    let input = "a / a";
    assert_eq!(format!("{}", expression(input)), "1")
}

#[cfg(test)]
fn expression(input: &str) -> Expr {
    let mut tokenizer = Tokenizer::new(input);
    let tokens = tokenizer.tokenize();
    let mut parser = Parser::new(tokens);
    parser.parse_expression().unwrap().simplify()
}