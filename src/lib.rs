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

pub fn evaluate_expr(expr: Expr) -> Expr {
    match &expr {
        Expr::Number(rational) => Expr::Number(rational.clone()),
        Expr::Symbol(s) => Expr::Symbol(s.clone()),
        Expr::Add(exprs) => {
            
            let mut symbol_count = 0;
            let mut sum = Rational::new(0, 1);
            let mut other_terms = Vec::new();
            for n in exprs.clone() {
                if let Expr::Symbol(_) = n {
                    symbol_count += 1;
                }
                if let Expr::Number(r) = n {
                    sum = sum + r;
                } else {
                    other_terms.push(evaluate_expr(n));
                }
            }
            
            let mut a = vec![Expr::Number(sum.clone())];
            a.append(&mut other_terms);

            

            if exprs.len() - 1 == symbol_count {
                if sum.numerator == 0 { 
                    return Expr::Add(other_terms);
                }
                return Expr::Add(a);
            }
            
            evaluate_expr(Expr::Add(a))

        },
        Expr::Sub(expr, expr1) => {
            
            match (*expr.clone(), *expr1.clone()) {
                (Expr::Number(n), Expr::Number(r)) => {
                    Expr::Number(n - r)
                },
                _=> {
                    Expr::Sub(Box::new(evaluate_expr(*expr.clone())), Box::new(evaluate_expr(*expr1.clone())))
                }
            }

        },
        Expr::Mul(exprs) => {
            let mut symbol_count = 0;
            let mut res = Rational::new(1, 1);
            let mut others = Vec::new();
            for n in exprs.clone() {
                if let Expr::Symbol(_) = n {
                    symbol_count += 1;
                }
                if let Expr::Number(m) = n {
                    res = res * m;
                } else {
                    others.push(n.clone());
                }
            }
            let mut res1 = vec![Expr::Number(res)];
            res1.append(&mut others);
            if exprs.len() - 1 == symbol_count {
                return Expr::Mul(res1);
            }
            evaluate_expr(Expr::Mul(res1))

        },
        Expr::Div(expr, expr1) => {
            match (*expr.clone(), *expr1.clone()) {
                (Expr::Number(n), Expr::Number(r)) => {
                    Expr::Number(n / r)
                },
                _=> {
                    Expr::Div(Box::new(evaluate_expr(*expr.clone())), Box::new(evaluate_expr(*expr1.clone())))
                }
            }
        },
        Expr::Pow(expr, expr1) => {
            todo!()
        },
        Expr::Neg(expr) => {
            match *expr.clone() {
                Expr::Number(n) => {
                    Expr::Number(-n)
                },
                _=> {
                    Expr::Neg(Box::new(evaluate_expr(*expr.clone())))
                }
            }
        }
    }
}
pub fn evaluate_expr1(expr: Expr) -> Expr {
    match expr {
        Expr::Number(rational) => Expr::Number(rational),
        Expr::Symbol(s) => Expr::Symbol(s),
        Expr::Add(exprs) => {
            let mut sum = Rational::new(0, 1);
            let mut symbols = HashMap::new();
            let mut other_terms = Vec::new();
            
            for expr in exprs {
                match evaluate_expr(expr) {
                    Expr::Number(r) => sum = sum + r,
                    Expr::Symbol(s) => {
                        *symbols.entry(s).or_insert(0) += 1;
                    },
                    Expr::Mul(factors) => {
                        // Handle coefficients like 2x, 3y, etc.
                        if let (coeff, symbol) = extract_coefficient(factors.clone()) {
                            if let Expr::Symbol(symbol) = symbol {
                                *symbols.entry(symbol).or_insert(0) += coeff.numerator;
                            }
                            
                        } else {
                            other_terms.push(Expr::Mul(factors));
                        }
                    },
                    other => other_terms.push(other),
                }
            }
            
            // Build the final expression
            let mut terms = Vec::new();
            
            // Add the constant term if non-zero
            if sum != Rational::new(0, 1) {
                terms.push(Expr::Number(sum));
            }
            
            // Add the symbols with their counts
            for (symbol, count) in symbols {
                match count {
                    1 => terms.push(Expr::Symbol(symbol)),
                    _ => terms.push(Expr::Mul(vec![
                        Expr::Number(Rational::new(count, 1)),
                        Expr::Symbol(symbol)
                    ])),
                }
            }
            
            // Add other terms
            terms.extend(other_terms);
            
            match terms.len() {
                0 => Expr::Number(Rational::new(0, 1)),
                1 => terms.into_iter().next().unwrap(),
                _ => Expr::Add(terms),
            }
        },
        Expr::Sub(a, b) => {
            let a_simpl = evaluate_expr(*a);
            let b_simpl = evaluate_expr(*b);
            
            match (a_simpl, b_simpl) {
                (Expr::Number(n1), Expr::Number(n2)) => Expr::Number(n1 - n2),
                (a, Expr::Number(n)) if n == Rational::new(0, 1) => a,
                (Expr::Add(terms), b) => {
                    // Distribute subtraction: (a + b + c) - d → a + b + c - d
                    let mut new_terms = terms;
                    new_terms.push(Expr::Neg(Box::new(b)));
                    evaluate_expr(Expr::Add(new_terms))
                },
                (a, b) => Expr::Sub(Box::new(a), Box::new(b)),
            }
        },
        Expr::Mul(exprs) => {
            let mut product = Rational::new(1, 1);
            let mut symbols = HashMap::new();
            let mut other_factors = Vec::new();
            
            for expr in exprs {
                match evaluate_expr(expr) {
                    Expr::Number(r) => product = product * r,
                    Expr::Symbol(s) => {
                        *symbols.entry(s).or_insert(0) += 1;
                    },
                    Expr::Pow(base, exp) => {
                        if let Expr::Symbol(s) = *base {
                            if let Expr::Number(n) = *exp {
                                *symbols.entry(s).or_insert(0) += n.numerator;
                            }
                            
                        } else {
                            other_factors.push(Expr::Pow(base, exp));
                        }
                    },
                    other => other_factors.push(other),
                }
            }
            
            // Build the final expression
            let mut factors = Vec::new();
            
            // Add the constant factor if not 1
            if product != Rational::new(1, 1) {
                factors.push(Expr::Number(product));
            }
            
            // Add the symbols with their exponents
            for (symbol, exp) in symbols {
                match exp {
                    1 => factors.push(Expr::Symbol(symbol)),
                    _ => factors.push(Expr::Pow(
                        Box::new(Expr::Symbol(symbol)),
                        Box::new(Expr::Number(Rational::new(exp, 1))),
                    )),
                }
            }
            
            // Add other factors
            factors.extend(other_factors);
            
            match factors.len() {
                0 => Expr::Number(Rational::new(1, 1)),
                1 => factors.into_iter().next().unwrap(),
                _ => Expr::Mul(factors),
            }
        },
        Expr::Div(a, b) => {
            let a_simpl = evaluate_expr(*a);
            let b_simpl = evaluate_expr(*b);
            
            match (a_simpl, b_simpl) {
                (Expr::Number(n), Expr::Number(m)) => Expr::Number(n / m),
                (a, Expr::Number(n)) if n == Rational::new(1, 1) => a,
                (a, b) => Expr::Div(Box::new(a), Box::new(b)),
            }
        },
        Expr::Pow(a, b) => {
            todo!();
            //let a_simpl = evaluate_expr(*a);
            //let b_simpl = evaluate_expr(*b);
            //
            //match (a_simpl, b_simpl) {
            //    (Expr::Number(n), Expr::Number(m)) => Expr::Number(n.pow(m.numerator())),
            //    (a, Expr::Number(n)) if n == Rational::new(1, 1) => a,
            //    (a, b) => Expr::Pow(Box::new(a), Box::new(b)),
            //}
        },
        Expr::Neg(a) => {
            match evaluate_expr(*a) {
                Expr::Number(n) => Expr::Number(-n),
                a_simpl => Expr::Neg(Box::new(a_simpl)),
            }
        },
        // Remove Expr::Res variant as it's not needed
    }
}


impl Expr {
    pub fn div(self, denominator: Expr) -> Expr {
        match (&self, &denominator) {
            (Expr::Symbol(n), Expr::Symbol(n1)) => {
                if n == n1 {
                    Expr::Number(Rational { numerator: 1, denominator: 1 })
                } else {
                    Expr::Div(Box::new(self), Box::new(denominator))
                }
            }

            _=> {Expr::Div(Box::new(self), Box::new(denominator))}
        }
        
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
                    
                    
                    Expr::Mul(v) => {
                        let mut simplified_factors = Vec::new();
                        let mut constant = Rational::new(1, 1);
                        
                        // First pass: simplify all factors and collect constants
                        for factor in v {
                            let simplified = factor.simplify();
                            match simplified {
                                Expr::Number(n) => constant = constant * n,
                                Expr::Mul(nested_factors) => {
                                    // Flatten nested multiplication
                                    for f in nested_factors {
                                        simplified_factors.push(f);
                                    }
                                }
                                _ => simplified_factors.push(simplified),
                            }
                        }
                        
                        // Add the combined constant if it's not 1
                        if constant != Rational::new(1, 1) {
                            simplified_factors.insert(0, Expr::Number(constant));
                        }
                        
                        // Second pass: combine like terms (x * x → x^2)
                        let mut combined_factors: HashMap<Expr, i64> = HashMap::new();
                        for factor in simplified_factors {
                            *combined_factors.entry(factor).or_insert(0) += 1;
                        }
                        
                        // Rebuild the factors
                        let mut new_factors = Vec::new();
                        for (base, count) in combined_factors {
                            if count == 1 {
                                new_factors.push(base);
                            } else {
                                new_factors.push(Expr::Pow(
                                    Box::new(base),
                                    Box::new(Expr::Number(Rational::new(count, 1)))
                                ));
                            }
                        }
                        
                        match new_factors.len() {
                            0 => Expr::Number(Rational::new(1, 1)), // Empty product is 1
                            1 => new_factors.into_iter().next().unwrap(),
                            _ => Expr::Mul(new_factors)
                        }
                    }
                    
                    Expr::Div(numerator, denominator) => {
                        let num = numerator.simplify();
                        let den = denominator.simplify();
        
                        // Case 1: a/a → 1
                        if num == den {
                            return Expr::Number(Rational::new(1, 1));
                        }
        
                        // Case 2: (a*b)/a → b
                        if let Expr::Mul(factors) = num.clone() {
                            if let Some(pos) = factors.iter().position(|f| f == &den) {
                                let mut new_factors = factors.clone();
                                new_factors.remove(pos);
                                return if new_factors.is_empty() {
                                    Expr::Number(Rational::new(1, 1))
                                } else if new_factors.len() == 1 {
                                    new_factors.into_iter().next().unwrap()
                                } else {
                                    Expr::Mul(new_factors)
                                };
                            }
                        }
        
                        // Case 3: (a/b)/c → a/(b*c)
                        if let Expr::Div(n, d) = num {
                            return Expr::Div(n, Box::new(Expr::Mul(vec![*d, den])));
                        }
        
                        // Default case
                        Expr::Div(Box::new(num), Box::new(den))
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
                    let needs_paren = false;
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