use std::collections::HashMap;

use crate::tokenizer::{Base, Rational};

pub mod parser;
pub mod tokenizer;
pub mod tests;


#[derive(Debug, PartialEq, Clone, Hash, Eq)]
pub enum Expr {
    Number(Rational),
    Symbol(String),
    Add(Vec<Expr>),
    Sub(Box<Expr>, Box<Expr>),
    Mul(Vec<Expr>),
    Div(Box<Expr>, Box<Expr>),
    Pow(Box<Expr>, Box<Expr>),
    Neg(Box<Expr>)
}

impl Expr {
    pub fn div(self, denominator: Expr) -> Expr {
        Expr::Div(Box::new(self), Box::new(denominator))
    }

    pub fn mult(self, factors: Vec<Expr>) -> Expr {
        let mut v = vec![self];
        v.append(&mut factors.clone());
        Expr::Mul(v)
    }

    pub fn simplify(self) -> Expr {
        match self {
            Expr::Add(terms) => {
                        // First simplify all sub-expressions
                        let mut simplified_terms: Vec<Expr> = terms.into_iter()
                            .map(|t| t.simplify())
                            .collect();
                
                        // Group like terms
                        let mut term_counts: HashMap<Expr, Rational> = HashMap::new();
                        let mut constant = Rational::new(0, 1);
                
                        for term in simplified_terms {
                            match term {
                                Expr::Number(n) => constant = constant+n,
                                Expr::Mul(factors) => {
                                    let (coeff, base) = extract_coefficient(factors);
                                    *term_counts.entry(base).or_default() = term_counts.get(&base).unwrap_or(&Rational::new(0,1)).clone() + coeff;
                                },
                                Expr::Symbol(s) => {
                                    let base = Expr::Symbol(s);
                                    *term_counts.entry(base).or_default() = term_counts.get(&base).unwrap_or(&Rational::new(0,1)).clone() + Rational::new(1,1);
                                },
                                _ => {
                                    // For other cases, just add to counts as-is
                                    *term_counts.entry(term).or_default() = term_counts.get(&term).unwrap_or(&Rational::new(0,1)).clone() + Rational::new(1,1);
                                }
                            }
                        }
                
                        // Rebuild the expression
                        let mut new_terms = Vec::new();
                
                        // Add the constant term if non-zero
                        if constant != Rational::new(0, 1) {
                            new_terms.push(Expr::Number(constant));
                        }
                
                        // Add the grouped variable terms
                        for (base, coeff) in term_counts {
                            if coeff == Rational::new(0, 1) {
                                continue;
                            }
                    
                            let term = if coeff == Rational::new(1, 1) {
                                base
                            } else {
                                Expr::Mul(vec![
                                    Expr::Number(coeff),
                                    base
                                ])
                            };
                            new_terms.push(term);
                        }
                
                        match new_terms.len() {
                            0 => Expr::Number(Rational::new(0, 1)),
                            1 => new_terms.into_iter().next().unwrap(),
                            _ => Expr::Add(new_terms)
                        }
                    },
                    Expr::Add(v) => {
                        let a = v[0].clone();
                        let b = v[1].clone();

                        let a = a.simplify();
                        let b = b.simplify();
                        match (a, b) {
                            (Expr::Number(x), Expr::Number(y)) => Expr::Number(x + y),
                            (Expr::Number(n), x) | (x, Expr::Number(n)) if n.numerator == 0 => x,
                            (a, b) => Expr::Add(vec![a, b]),
                        }
                    }
                    
                    Expr::Sub(a, b) => {
                        Expr::Add(vec![a.simplify(), Expr::Neg(Box::new(b.simplify())).simplify()])
                    }
                    
                    Expr::Mul(v) => {
                        let a = v[0].clone();
                        if v.get(1).is_none() {
                            return a.simplify();
                        }
                        let b = v[1].clone();
                        let a = a.simplify();
                        let b = b.simplify();
                        match (a, b) {
                            (Expr::Number(x), Expr::Number(y)) => Expr::Number(x * y),
                            (Expr::Number(n), x) | (x, Expr::Number(n)) => match n.numerator {
                                0 => Expr::Number(Rational::new(0, 1)),
                                1 => x,
                                _ => Expr::Mul(vec![Expr::Number(n), x]),
                            },
                            (a, b) => Expr::Mul(vec![a, b]),
                        }
                    }
                    
                    Expr::Div(a, b) => {
                        let a = a.simplify();
                        let b = b.simplify();
                        
                        match (a, b) {
                            // Number division
                            (Expr::Number(num), Expr::Number(den)) => {
                                Expr::Number(Rational::new(
                                    num.numerator * den.denominator,
                                    num.denominator * den.numerator,
                                ))
                            },
                            
                            // x / 1 → x
                            (expr, Expr::Number(den)) if den == Rational::new(1, 1) => expr,
                            
                            // (a * b) / c → a * (b / c)
                            (Expr::Mul(factors), denom) => {
                                let last = factors.last().unwrap().clone();
                                let rest = factors[..factors.len()-1].to_vec();
                                Expr::Mul(rest).mult(vec![last.div(denom)])
                            },
                            
                            // Default case
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
                    Expr::Add(v) if v[0] == Expr::Number(Rational::new(0, 1)) => v[1].clone().simplify(),
                    Expr::Add(v) if v[1] == Expr::Number(Rational::new(0, 1)) => v[0].clone().simplify(),
                    Expr::Mul(v) if v[0] == Expr::Number(Rational::new(1, 1)) => v[1].clone().simplify(),
                    Expr::Mul(v) if v[1] == Expr::Number(Rational::new(1, 1)) => v[0].clone().simplify(),
                    
                    // Recursively simplify other cases
                    other => other,
        }
    }
}

fn extract_coefficient(factors: Vec<Expr>) -> (Rational, Expr) {
    let mut coeff = Rational::new(1, 1);
    let mut other_factors = Vec::new();
    
    for factor in factors {
        match factor {
            Expr::Number(n) => coeff = coeff * n,
            _ => other_factors.push(factor)
        }
    }
    
    let base = if other_factors.is_empty() {
        Expr::Number(Rational::new(1, 1))
    } else if other_factors.len() == 1 {
        other_factors.into_iter().next().unwrap()
    } else {
        Expr::Mul(other_factors)
    };
    
    (coeff, base)
}

fn flatten_expr(expr: Expr) -> Expr {
    match expr {
        Expr::Add(v) => {
            let a = v.first().unwrap();
            let b = &v[1];
            let mut terms = Vec::new();
            terms.extend(flatten_add(a.clone()));
            terms.extend(flatten_add(b.clone()));
            Expr::Add(terms)
        },
        Expr::Mul(v) => {
            let a = v.first().unwrap();
            let b = &v[1];
            let mut factors = Vec::new();
            factors.extend(flatten_mul(a.clone()));
            factors.extend(flatten_mul(b.clone()));
            Expr::Mul(factors)
        },
        other => other
    }
}

fn flatten_add(expr: Expr) -> Vec<Expr> {
    match flatten_expr(expr) {
        Expr::Add(terms) => terms,
        other => vec![other]
    }
}

fn flatten_mul(expr: Expr) -> Vec<Expr> {
    match flatten_expr(expr) {
        Expr::Mul(factors) => factors,
        other => vec![other]
    }
}

pub fn add(left: u64, right: u64) -> u64 {
    left + right
}


use std::fmt;

impl fmt::Display for Expr {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Expr::Number(r) => write!(f, "{}", r),
            Expr::Symbol(s) => write!(f, "{}", s),
            
            Expr::Add(terms) => {
                let mut first = true;
                for term in terms {
                    let needs_paren = matches!(term, Expr::Sub(..) | Expr::Neg(..));
                    write!(f, "{}{}", 
                        if first { "" } else { " + " }, 
                        if needs_paren { format!("({})", term) } else { format!("{}", term) }
                    )?;
                    first = false;
                }
                Ok(())
            },
            
            Expr::Neg(expr) => write!(f, "-{}", 
                if expr.needs_paren() { format!("({})", expr) } else { format!("{}", expr) }
            ),
            
            Expr::Sub(a, b) => write!(f, "{} - {}", 
                if a.needs_paren() { format!("({})", a) } else { format!("{}", a) },
                if b.needs_paren() { format!("({})", b) } else { format!("{}", b) }
            ),
            
            Expr::Mul(factors) => {
                let mut first = true;
                for factor in factors {
                    write!(f, "{}{}", 
                        if first { "" } else { " * " },
                        if factor.needs_paren_mul() { format!("({})", factor) } else { format!("{}", factor) }
                    )?;
                    first = false;
                }
                Ok(())
            },
            
            Expr::Div(a, b) => write!(f, "{}/{}", 
                if a.needs_paren_div() { format!("({})", a) } else { format!("{}", a) },
                if b.needs_paren_div() { format!("({})", b) } else { format!("{}", b) }
            ),
            
            Expr::Pow(a, b) => write!(f, "{}^{}", 
                if a.needs_paren_pow() { format!("({})", a) } else { format!("{}", a) },
                if b.needs_paren_pow() { format!("({})", b) } else { format!("{}", b) }
            ),
        }
    }
}

// Helper trait for parenthesis decisions
trait NeedsParen {
    fn needs_paren(&self) -> bool;
    fn needs_paren_mul(&self) -> bool;
    fn needs_paren_div(&self) -> bool;
    fn needs_paren_pow(&self) -> bool;
}

impl NeedsParen for Expr {
    fn needs_paren(&self) -> bool {
        matches!(self, Expr::Add(_) | Expr::Sub(..) | Expr::Neg(_))
    }
    
    fn needs_paren_mul(&self) -> bool {
        matches!(self, Expr::Add(_) | Expr::Sub(..) | Expr::Div(..))
    }
    
    fn needs_paren_div(&self) -> bool {
        matches!(self, Expr::Add(_) | Expr::Sub(..) | Expr::Mul(_) | Expr::Div(..))
    }
    
    fn needs_paren_pow(&self) -> bool {
        matches!(self, Expr::Add(_) | Expr::Sub(..) | Expr::Mul(_) | Expr::Div(..) | Expr::Pow(..))
    }
}