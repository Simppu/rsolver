use crate::tokenizer::{Base, Rational};

pub mod parser;
pub mod tokenizer;
pub mod tests;


#[derive(Debug, PartialEq, Clone)]
pub enum Expr {
    Number(Rational),
    Symbol(String),
    Add(Box<Expr>, Box<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Box<Expr>, Box<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>)
}

impl Expr {
    fn contains_symbol(&self, sym: &str) -> bool {
        match self {
            Expr::Symbol(s) => s == sym,
            Expr::Add(a, b) | Expr::Sub(a, b) | Expr::Mul(a, b) | Expr::Div(a, b) => {
                a.contains_symbol(sym) || b.contains_symbol(sym)
            }
            Expr::Neg(e) => e.contains_symbol(sym),
            _ => false,
        }
    }
    pub fn simplify(self) -> Expr {
        match self {
            

            // Basic arithmetic simplification
            Expr::Add(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (a.clone(), b) {
                    (Expr::Number(x), Expr::Number(y)) => Expr::Number(x + y),
                    (Expr::Symbol(ref s1), Expr::Symbol(ref s2)) if s1 == s2 => {
                        Expr::Mul(
                            Box::new(Expr::Number(Rational::new(2, 1))),
                            Box::new(Expr::Symbol(s1.clone()))
                        )
                    },
                    (Expr::Number(n), x) | (x, Expr::Number(n)) if n.numerator == 0 => x,
                    (Expr::Mul(c1, s1), Expr::Mul(c2, s2))
                     if *s1 == *s2 => {
                        Expr::Mul(
                            Box::new(Expr::Add(c1, c2).simplify()),
                            Box::new(*s1.clone())
                        )
                    },
                    
                    // a + a*x → (1 + x)*a when x is complex
                    (Expr::Symbol(ref s), Expr::Mul(c, b)) | 
                    (Expr::Mul(c, b), Expr::Symbol(ref s))
                    if b.contains_symbol(s) => {
                        // Advanced grouping would go here
                        Expr::Add(Box::new(a), Box::new(*b))
                    },
                    
                    // Default case
                    (a, b) => Expr::Add(Box::new(a), Box::new(b)),
                    
                }
            }
            
            Expr::Sub(a, b) => {
                Expr::Add(Box::new(a.simplify()), Box::new(Expr::Neg(Box::new(b.simplify())))).simplify()
            }
            
            Expr::Mul(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (a, b) {
                    (Expr::Symbol(ref s1), Expr::Symbol(ref s2)) if s1 == s2 => {
                        Expr::Pow(
                            Box::new(Expr::Symbol(s1.clone())),
                            Box::new(Expr::Number(Rational::new(2, 1)))
                        )
                    },
                    (Expr::Number(x), Expr::Number(y)) => Expr::Number(x * y),
                    (Expr::Number(n), x) | (x, Expr::Number(n)) => match n.numerator {
                        0 => Expr::Number(Rational::new(0, 1)),
                        1 => x,
                        _ => Expr::Mul(Box::new(Expr::Number(n)), Box::new(x)),
                    },
                    (a, b) => Expr::Mul(Box::new(a), Box::new(b)),
                }
            }
            
            Expr::Div(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (a, b) {
                    (Expr::Number(x), Expr::Number(y)) => {
                        Expr::Number(Rational::new(
                            x.numerator * y.denominator,
                            x.denominator * y.numerator,
                        ))
                    }
                    (a, Expr::Number(n)) if n.numerator == 1 => {
                        Expr::Mul(Box::new(a), Box::new(Expr::Number(Rational::new(n.denominator, 1))))
                    }
                    (a, b) => Expr::Div(Box::new(a), Box::new(b)),
                }
            }
            
            Expr::Neg(a) => {
                let a = a.simplify();
                match a {
                    Expr::Number(n) => Expr::Number(-n),
                    Expr::Neg(x) => *x, // Double negative
                    a => Expr::Neg(Box::new(a)),
                }
            }
            
            Expr::Pow(a, b) => {
                let a = a.simplify();
                let b = b.simplify();
                match (a, b) {
                    (Expr::Number(x), Expr::Number(y)) if y.denominator == 1 => {
                        Expr::Number(Rational::new(
                            x.numerator.pow(y.numerator as u32),
                            x.denominator.pow(y.numerator as u32),
                        ))
                    }
                    (a, Expr::Number(n)) if n.numerator == 1 && n.denominator == 1 => a,
                    (a, b) => Expr::Pow(Box::new(a), Box::new(b)),
                }
            }
            
            // Symbolic simplification rules
            Expr::Add(a, b) if *a == Expr::Number(Rational::new(0, 1)) => b.simplify(),
            Expr::Add(a, b) if *b == Expr::Number(Rational::new(0, 1)) => a.simplify(),
            Expr::Mul(a, b) if *a == Expr::Number(Rational::new(1, 1)) => b.simplify(),
            Expr::Mul(a, b) if *b == Expr::Number(Rational::new(1, 1)) => a.simplify(),
            
            // Recursively simplify other cases
            other => other,
        }
    }
}



pub fn add(left: u64, right: u64) -> u64 {
    left + right
}


