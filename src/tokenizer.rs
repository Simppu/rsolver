use std::{fmt::Display, ops::{Add, Div, Mul, Neg, Sub}, str::Chars};


#[derive(Debug, PartialEq, Clone)]
pub enum Token {
    Number(Rational),
    Symbol(String),
    Plus,
    Minus,
    Star,
    Slash,
    Caret,
    LParen,
    RParen,
    Eq,
    Log(Base),
    Eof,
}

#[derive(Debug, PartialEq, Clone)]
pub enum Base {
    N,
    B,
    Lg,
    Custom(i64)
}

#[allow(clippy::derived_hash_with_manual_eq)]
#[derive(Debug, Clone, Hash, Eq)]
pub struct Rational {
    pub numerator: i64,
    pub denominator: i64, // Always positive
}

impl Default for Rational {
    fn default() -> Self {
        Self { numerator: 0, denominator: 1 }
    }
}

impl Display for Rational {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.denominator == 1 {
            write!(f, "{}", self.numerator)
        } else {
            write!(f, "{}/{}", self.numerator, self.denominator)
        }
    }
}


impl Rational {
    pub fn new(n: i64, d: i64) -> Self {
        if d == 1 {
            return Self {
                numerator: n,
                denominator: 1
            };
        }
        let gcd = num::integer::gcd(n.abs(), d.abs());
        let sign = if d < 0 { -1 } else { 1 };
        Rational {
            numerator: sign * n / gcd,
            denominator: (d * sign) / gcd,
        }
    }

    fn invert(&self) -> Self {
        Self { numerator: self.denominator, denominator: self.numerator }
    }
}

impl Add for Rational {
    type Output = Rational;

    fn add(self, rhs: Self) -> Self::Output {
        let n = self.numerator * rhs.denominator + rhs.numerator * self.denominator;
        let d = self.denominator * rhs.denominator;

        let gcd = num::integer::gcd(n.abs(), d.abs());
        let sign = if d < 0 { -1 } else {1};

        Self {
            numerator: sign * n / gcd,
            denominator: (d * sign) / gcd
        }
    }
}

impl Neg for Rational {
    type Output = Rational;

    fn neg(self) -> Self::Output {
        Self {
            numerator: -self.numerator,
            denominator: self.denominator
        }
    }
}

impl PartialEq for Rational {
    fn eq(&self, other: &Self) -> bool {
        self.numerator * other.denominator == self.denominator * other.numerator
    }
}

impl PartialOrd for Rational {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let n1 = self.numerator * other.denominator;
        let n2 = other.numerator * self.denominator;
        
        n1.partial_cmp(&n2)
        
    }
}

impl Sub for Rational {
    type Output = Rational;

    fn sub(self, rhs: Self) -> Self::Output {
        let n = self.numerator * rhs.denominator - rhs.numerator * self.denominator;
        let d = self.denominator * rhs.denominator;

        let gcd = num::integer::gcd(n.abs(), d.abs());
        let sign = if d < 0 { -1 } else {1};

        Self {
            numerator: sign * n / gcd,
            denominator: (d * sign) / gcd
        }
    }
}

impl Mul for Rational {
    type Output = Rational;

    fn mul(self, rhs: Self) -> Self::Output {
        let n = self.numerator * rhs.numerator;
        let d = self.denominator * rhs.denominator;

        let gcd = num::integer::gcd(n.abs(), d.abs());
        let sign = if d < 0 { -1 } else {1};

        Self {
            numerator: sign * n / gcd,
            denominator: (d * sign) / gcd
        }
    }
}

impl Div for Rational {
    type Output = Rational;

    fn div(self, rhs: Self) -> Self::Output {
        self * rhs.invert()
    }
}


pub struct Tokenizer<'a> {
    chars: Chars<'a>,
    look_ahead: Option<char>
}

impl<'a> Tokenizer<'a> {
    pub fn new(input: &'a str) -> Self {
        let mut chars = input.chars();
        let look_ahead = chars.next();

        Tokenizer { chars, look_ahead }
    }

    fn next_char(&mut self) -> Option<char> {
        let current = self.look_ahead;
        self.look_ahead = self.chars.next();
        current
    }

    fn peek(&self) -> Option<char> {
        self.look_ahead
    }

    fn skip_whitespace(&mut self) {
        while let Some(c) = self.peek() {
            if !c.is_whitespace() {
                break;
            }
            self.next_char();
        }
    }

    


    fn read_number(&mut self, first: char) -> Result<Token, String> {
        let mut s = String::new();
        s.push(first);
        while let Some(c) = self.peek() {
            match c {
                '0'..='9' | '/' => s.push(self.next_char().unwrap()),
                _ => break,
            }
        }
        


        if let Some(pos) = s.find('/') {
            let num = s[..pos].parse::<i64>().map_err(|e| e.to_string())?;
            let den = s[pos+1..].parse::<i64>().map_err(|e| e.to_string())?;
            Ok(Token::Number(Rational::new(num, den)))
        } else {
            let num = s.parse::<i64>().map_err(|e| e.to_string())?;
            Ok(Token::Number(Rational::new(num, 1)))
        }
    }


    fn read_symbol(&mut self, first: char) -> Token {
        let mut s = String::new();
        s.push(first);

        while let Some(c) = self.peek() {
            if c.is_alphabetic() {
                s.push(self.next_char().unwrap());
            } else {
                break;
            }
        }

        match s.as_str() {
            "log" => {
                Token::Log(Base::Lg)
            },
            "ln" => {
                Token::Log(Base::N)
            },
            "lb" => {Token::Log(Base::B)},
            _=> {
                Token::Symbol(s)
            }
        }

        
    }

    pub fn next_token(&mut self) -> Token {
        self.skip_whitespace();

        match self.next_char() {
            Some('+') => Token::Plus,
            Some('-') => Token::Minus,
            Some('*') => Token::Star,
            Some('/') => Token::Slash,
            Some('^') => Token::Caret,
            Some('(') => Token::LParen,
            Some(')') => Token::RParen,
            Some('=') => Token::Eq,
            Some(c) if c.is_ascii_digit() => self.read_number(c).unwrap(),
            Some(c) if c.is_alphabetic() => self.read_symbol(c),
            None => Token::Eof,
            Some(c) => panic!("Unexpected character: {}", c),
        }
    }

    pub fn tokenize(&mut self) -> Vec<Token> {
        let mut tokens = Vec::new();
        loop {
            let token = self.next_token();
            if token == Token::Eof {break;}
            tokens.push(token);
        }

        tokens
    }


}
