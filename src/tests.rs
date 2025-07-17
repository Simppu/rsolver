use crate::tokenizer::{Rational, Token, Tokenizer};



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