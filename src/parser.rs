use std::{collections::HashMap, iter::Peekable, string::ParseError};

use crate::{tokenizer::Token, Expr};




fn parse_expr(input: &str) -> Result<Expr, ParseError> {
    todo!()
}







//pub fn evaluate(node: &Expr, vars: &HashMap<String, f64>) -> f64 {
//    match node {
//        Expr::Number(n) => *n,
//        Expr::Symbol(s) => vars[s],
//        Expr::Add(a, b) => evaluate(a, vars) + evaluate(b, vars),
//        Expr::Sub(a, b) => evaluate(a, vars) - evaluate(b, vars),
//        Expr::Mul(a, b) => evaluate(a, vars) * evaluate(b, vars),
//        Expr::Div(a, b) => evaluate(a, vars) / evaluate(b, vars),
//        Expr::Pow(a, b) => evaluate(a, vars).powf(evaluate(b, vars)),
//        Expr::Neg(a) => - evaluate(a, vars),
//    }
//}



pub struct Parser {
    tokens: Peekable<std::vec::IntoIter<Token>>,
    current: Option<Token>
}

impl Parser{
    pub fn new(tokens: Vec<Token>) -> Self {
        let mut peekable = tokens.into_iter().peekable();
        let current = peekable.next();
        Parser {
            tokens: peekable,
            current,
        }
    }

    fn advance(&mut self) {
        self.current = self.tokens.next();
    }

    fn peek(&mut self) -> Option<&Token> {
        self.tokens.peek()
    }

    fn consume(&mut self, expected: Token) -> Result<(), String> {
        match &self.current {
            Some(token) if token == &expected => {
                self.advance();
                Ok(())
            }
            Some(other) => Err(format!(
                "Expected {:?}, found {:?}",
                expected, other
            )),
            None => Err(format!("Expected {:?}, found EOF", expected)),
        }
    }
}


impl Parser{
    pub fn parse_expression(&mut self) -> Result<Expr, String> {
        self.parse_additive()
    }

    fn parse_additive(&mut self) -> Result<Expr, String> {
        let mut left = self.parse_multiplicative()?;
        
        while let Some(token) = &self.current {
            match token {
                Token::Plus | Token::Minus => {
                    let op = token.clone();
                    self.advance();
                    let right = self.parse_multiplicative()?;
                    
                    left = match op {
                        Token::Plus => Expr::Add(vec![left, right]),
                        Token::Minus => Expr::Sub(Box::new(left), Box::new(right)),
                        _ => unreachable!(),
                    };
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_multiplicative(&mut self) -> Result<Expr, String> {
        let mut left = self.parse_exponent()?;
        
        while let Some(token) = &self.current {
            match token {
                Token::Star => {
                    self.advance();
                    let right = self.parse_exponent()?;
                    left = Expr::Mul(vec![left, right]);
                }
                Token::Slash => {
                    self.advance();
                    let right = self.parse_exponent()?;
                    left = Expr::Div(Box::new(left), Box::new(right));
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_exponent(&mut self) -> Result<Expr, String> {
        let mut left = self.parse_primary()?;
        
        while let Some(Token::Caret) = &self.current {
            self.advance();
            let right = self.parse_primary()?;
            left = Expr::Pow(Box::new(left), Box::new(right));
        }
        Ok(left)
    }

    fn parse_primary(&mut self) -> Result<Expr, String> {
        match self.current.take() {
            Some(Token::Number(n)) => {
                self.advance();
                Ok(Expr::Number(n))
            }
            Some(Token::Symbol(s)) => {
                self.advance();
                Ok(Expr::Symbol(s))
            }
            Some(Token::LParen) => {
                self.advance();
                let expr = self.parse_expression()?;
                self.consume(Token::RParen)?;
                Ok(expr)
            }
            Some(Token::Minus) => {
                self.advance();
                Ok(Expr::Neg(Box::new(self.parse_primary()?)))
            }
            Some(t) => Err(format!("Unexpected token: {:?}", t)),
            None => Err("Unexpected end of input".into()),
        }
    }
}