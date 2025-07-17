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
        let terms = self.parse_additive()?;
        Ok(if terms.len() == 1 {
            terms.into_iter().next().unwrap()
        } else {
            Expr::Add(terms)
        })
    }

    fn parse_additive(&mut self) -> Result<Vec<Expr>, String> {
        let mut left = self.parse_multiplicative()?;
        
        while let Some(token) = &self.current {
            match token {
                Token::Plus | Token::Minus => {
                    let op = token.clone();
                    self.advance();
                    let mut right = self.parse_multiplicative()?;
                    
                    if op == Token::Minus {
                        right = right.into_iter()
                            .map(|expr| Expr::Neg(Box::new(expr)))
                            .collect();
                    }
                    
                    left.extend(right);
                }
                _ => break,
            }
        }
        Ok(left)
    }

    fn parse_multiplicative(&mut self) -> Result<Vec<Expr>, String> {
        let mut left = self.parse_exponent()?;
        
        while let Some(token) = &self.current {
            match token {
                Token::Star | Token::Slash => {
                    let op = token.clone();
                    self.advance();
                    let mut right = self.parse_exponent()?;
                    
                    if op == Token::Slash {
                        right = vec![Expr::Div(
                            Box::new(Expr::Mul(right)),
                            Box::new(Expr::Mul(left)),
                        )];
                        left = right;
                        break;
                    }
                    
                    left.extend(right);
                }
                _ => break,
            }
        }
        Ok(left)
    }


    fn parse_exponent(&mut self) -> Result<Vec<Expr>, String> {
        let mut factors = self.parse_primary()?;
        
        while let Some(Token::Caret) = &self.current {
            self.advance();
            let exponent = self.parse_primary()?;
            factors = vec![Expr::Pow(
                Box::new(Expr::Mul(factors)),
                Box::new(Expr::Mul(exponent)),
            )];
        }
        Ok(factors)
    }

    fn parse_primary(&mut self) -> Result<Vec<Expr>, String> {
        match self.current.take() {
            Some(Token::Number(n)) => {
                self.advance();
                Ok(vec![Expr::Number(n)])
            }
            Some(Token::Symbol(s)) => {
                self.advance();
                Ok(vec![Expr::Symbol(s)])
            }
            Some(Token::LParen) => {
                self.advance();
                let expr = self.parse_expression()?;
                self.consume(Token::RParen)?;
                Ok(vec![expr])
            }
            Some(Token::Minus) => {
                self.advance();
                Ok(self.parse_primary()?
                    .into_iter()
                    .map(|expr| Expr::Neg(Box::new(expr)))
                    .collect())
            }
            Some(t) => Err(format!("Unexpected token: {:?}", t)),
            None => Err("Unexpected end of input".into()),
        }
    }
}