// Since we might need to parse large trees, we'd prefer to write a lexer that
// almost never allocates Strings.  The easy way to do this is to force the
// entire input to be an in-memory string slice (`&str`) and then tie the lifetime
// of any returned string slices to that of the input.  This is the solution used by the
// [Rust compiler's lexer](https://doc.rust-lang.org/1.90.0/nightly-rustc/rustc_lexer/index.html),
// for example.  We'll do that for now, but eventually, I think it might be possible to use
// a `String` that gets cleared and reused.  The trick is to pass that `String` by value to
// each call to the lexer, then the lexer returns that `String` again by value to the caller.
// That way, the caller gives up ownership of the `String` at each call to the lexer,
// which can only happen if there are no active borrows at the time.  Food for thought...

// TODO:
// - Poison lexer once an error has been returned
// - Add position indicators to errors and remove `char` from LexerError
// - Add documentation comments

use std::collections::HashMap;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Token<'a> {
    // Chunky tokens
    Attributes(&'a str),
    Number(f64, &'a str), // Preserve original characters in case used as label, e.g., `00024`
    UnquotedLabel(&'a str),
    QuotedLabel(&'a str),

    // Symbols
    LeftParen,  // (
    Comma,      // ,
    RightParen, // )
    Colon,      // :
    Semicolon,  // ;

    // Marker
    Eof, // End of file
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum LexerError {
    InvalidCharacter(char),
    InvalidNumber,
    UnterminatedComment,
    UnterminatedString,
}

pub struct Lexer<'a> {
    input: &'a str,
    pos: usize,
}

impl<'a> Lexer<'a> {
    pub fn new(input: &'a str) -> Self {
        Lexer { input, pos: 0 }
    }

    fn peek(&self) -> Option<char> {
        self.input[self.pos..].chars().next()
    }

    fn advance(&mut self) -> Option<char> {
        match self.peek() {
            None => None,
            Some(c) => {
                self.pos += c.len_utf8();
                Some(c)
            }
        }
    }

    fn consume_and<T>(&mut self, t: T) -> T {
        self.advance();
        t
    }

    pub fn next_token(&mut self) -> Result<Token<'a>, LexerError> {
        self.skip_whitespace();
        match self.peek() {
            None => Ok(Token::Eof),
            Some(c) => match c {
                '(' => self.consume_and(Ok(Token::LeftParen)),
                ')' => self.consume_and(Ok(Token::RightParen)),
                ',' => self.consume_and(Ok(Token::Comma)),
                ':' => self.consume_and(Ok(Token::Colon)),
                ';' => self.consume_and(Ok(Token::Semicolon)),

                ']' => self.consume_and(Err(LexerError::InvalidCharacter(c))),
                '[' => {
                    // ATTRIBUTES or COMMENT
                    match self.lex_comment_or_attribute() {
                        Ok(contents) => {
                            if contents.starts_with("[&") {
                                Ok(Token::Attributes(contents))
                            } else {
                                self.next_token() // Comment
                            }
                        }
                        Err(error) => Err(error),
                    }
                }

                '\'' => self.lex_quoted_label(),

                c => {
                    if c.is_ascii_digit() || c == '-' || c == '+' {
                        self.lex_number()
                    } else {
                        self.lex_unquoted_label()
                    }
                }
            },
        }
    }

    fn maybe_match_if(&mut self, predicate: impl FnOnce(char) -> bool) -> bool {
        match self.peek() {
            Some(c) if predicate(c) => {
                self.advance();
                true
            }
            _ => false,
        }
    }

    fn maybe_match_char(&mut self, target: char) -> bool {
        self.maybe_match_if(|c| c == target)
    }

    fn match_digits_plus(&mut self) -> bool {
        if !self.maybe_match_if(|c| c.is_ascii_digit()) {
            false
        } else {
            while self.maybe_match_if(|c| c.is_ascii_digit()) {}
            true
        }
    }

    //WS             : [ \r\n\t]+ -> skip ;
    fn skip_whitespace(&mut self) {
        while self.maybe_match_if(|c| c.is_ascii_whitespace()) {}
    }

    //NUMBER         : [+-]? [0-9]+ ('.' [0-9]+)? ([eE] [+-]? [0-9]+)? ;
    fn lex_number(&mut self) -> Result<Token<'a>, LexerError> {
        let start = self.pos;

        self.maybe_match_if(|c| c == '+' || c == '-');
        if !self.match_digits_plus() {
            return Err(LexerError::InvalidNumber);
        }
        if self.maybe_match_char('.') {
            if !self.match_digits_plus() {
                return Err(LexerError::InvalidNumber);
            }
        }
        if self.maybe_match_if(|c| c == 'e' || c == 'E') {
            self.maybe_match_if(|c| c == '+' || c == '-');
            if !self.match_digits_plus() {
                return Err(LexerError::InvalidNumber);
            }
        }

        let s = &self.input[start..self.pos];
        let num = s.parse().map_err(|_| LexerError::InvalidNumber)?;
        Ok(Token::Number(num, s))
    }

    // // ATTRIBUTES comes before COMMENT so that "[&hello]" is recognized
    // // as an ATTRIBUTES instead of a COMMENT
    // ATTRIBUTES     : '[&' .*? ']' ;
    // COMMENT        : '[' .*? ']' -> skip ;
    fn lex_comment_or_attribute(&mut self) -> Result<&'a str, LexerError> {
        let start = self.pos;
        assert_eq!(
            self.advance(),
            Some('['),
            "lex_comment_or_attribute only called if lookahead matches ["
        );
        loop {
            match self.advance() {
                Some(']') => break,
                Some(_) => (),
                None => return Err(LexerError::UnterminatedComment),
            }
        }
        Ok(&self.input[start..self.pos])
    }

    fn lex_quoted_label(&mut self) -> Result<Token<'a>, LexerError> {
        let start = self.pos;
        assert_eq!(
            self.advance(),
            Some('\''),
            "lex_quoted_label only called if lookahead matches \'"
        );
        loop {
            match self.advance() {
                None => return Err(LexerError::UnterminatedString),
                Some('\'') => {
                    if !self.maybe_match_char('\'') {
                        break;
                    }
                }
                _ => (),
            }
        }
        Ok(Token::QuotedLabel(&self.input[start..self.pos]))
    }

    fn lex_unquoted_label(&mut self) -> Result<Token<'a>, LexerError> {
        let start = self.pos;
        loop {
            match self.peek() {
                None => break,
                Some(c) => match c {
                    ' ' | '\r' | '\n' | '\t' | '(' | ')' | '[' | ']' | '\'' | ':' | ';' | ',' => {
                        break;
                    }
                    _ => self.advance(),
                },
            };
        }

        if start == self.pos {
            let c = self.peek().expect("lex_unquoted_label should only have been called when there was a lookahead character");
            Err(LexerError::InvalidCharacter(c))
        } else {
            Ok(Token::UnquotedLabel(&self.input[start..self.pos]))
        }
    }
}

impl<'a> Iterator for Lexer<'a> {
    type Item = Result<Token<'a>, LexerError>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.next_token() {
            Ok(Token::Eof) => None,
            tok @ _ => Some(tok),
        }
    }
}

pub fn unquote(s: &str) -> Option<String> {
    if s.len() < 2 {
        None
    } else {
        let mut chars = s.chars();
        // Skip leading quote
        if chars.next() != Some('\'') {
            return None;
        }
        // Skip trailing quote
        if chars.next_back() != Some('\'') {
            return None;
        }
        let mut chars = chars.peekable();

        let mut result = String::new();
        while let Some(c) = chars.next() {
            result.push(c);
            if c == '\'' && chars.peek() == Some(&'\'') {
                // Skip quoted single quote
                chars.next();
            }
        }
        Some(result)
    }
}

pub fn parse_attributes<'a>(s: &'a str, attrs: &mut HashMap<String, String>) {
    let mut chars = s.chars();
    assert_eq!(chars.next(), Some('['));
    assert_eq!(chars.next(), Some('&'));
    assert_eq!(chars.next_back(), Some(']'));
    let mut done = false;
    while !done {
        // Every iteration of this loop consists of one new attribute
        let mut key = String::new();
        let mut value = String::new();
        let mut has_value = false;

        loop {
            match chars.next() {
                None => {
                    done = true;
                    break;
                }
                Some('=') => {
                    has_value = true;
                    break;
                }
                Some(';') => break,
                Some(c) => key.push(c),
            }
        }
        if key == "" {
            break;
        }
        if has_value {
            loop {
                match chars.next() {
                    None => {
                        done = true;
                        break;
                    }
                    Some(';') => break,
                    Some(c) => value.push(c),
                }
            }
        }

        attrs.insert(key, value);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {
        let test_str = "[Hello] [&R](A:0.1,'E(''2'')':[&mutations={A1C, 0.2}]0.5)";
        let tokens: Vec<_> = Lexer::new(test_str).collect();
        assert_eq!(
            tokens,
            vec![
                Ok(Token::Attributes("[&R]")),
                Ok(Token::LeftParen),
                Ok(Token::UnquotedLabel("A")),
                Ok(Token::Colon),
                Ok(Token::Number(0.1, "0.1")),
                Ok(Token::Comma),
                Ok(Token::QuotedLabel("'E(''2'')'")),
                Ok(Token::Colon),
                Ok(Token::Attributes("[&mutations={A1C, 0.2}]")),
                Ok(Token::Number(0.5, "0.5")),
                Ok(Token::RightParen),
            ]
        );
    }

    #[test]
    fn unquote_test() {
        // Ok cases
        assert_eq!(unquote("''"), Some(String::from("")));
        assert_eq!(unquote("'Hi'"), Some(String::from("Hi")));
        assert_eq!(unquote("'McDonald''s'"), Some(String::from("McDonald's")));

        // Bad cases
        assert_eq!(unquote(""), None);
        assert_eq!(unquote("A"), None);
        assert_eq!(unquote("Unquoted"), None);
        assert_eq!(unquote("'A"), None);
        assert_eq!(unquote("A'"), None);

        // Borderline cases that we recover from as best we can (lexer should never recognize these)
        //
        // dangling trailing single quote
        assert_eq!(unquote("'A''"), Some(String::from("A'")));
        // unquoted single quote
        assert_eq!(unquote("'McDonald's'"), Some(String::from("McDonald's")));
    }

    #[test]
    fn parse_attributes_empty() {
        let mut result = HashMap::new();
        parse_attributes("[&]", &mut result);
        assert_eq!(result, HashMap::new());
    }

    #[test]
    fn parse_attributes_one_single() {
        let mut result = HashMap::new();
        parse_attributes("[&apple]", &mut result);
        assert_eq!(
            result,
            HashMap::from([(String::from("apple"), String::from("")),])
        );
    }

    #[test]
    fn parse_attributes_two_single() {
        let mut result = HashMap::new();
        parse_attributes("[&apple;banana]", &mut result);
        assert_eq!(
            result,
            HashMap::from([
                (String::from("apple"), String::from("")),
                (String::from("banana"), String::from("")),
            ])
        );
    }

    #[test]
    fn parse_attributes_one_kv() {
        let mut result = HashMap::new();
        parse_attributes("[&apple=orange]", &mut result);
        assert_eq!(
            result,
            HashMap::from([(String::from("apple"), String::from("orange")),])
        );
    }

    #[test]
    fn parse_attributes_two_kv() {
        let mut result = HashMap::new();
        parse_attributes("[&apple=orange;mutations={A23C,5.6,C1T,1.2}]", &mut result);
        assert_eq!(
            result,
            HashMap::from([
                (String::from("apple"), String::from("orange")),
                (String::from("mutations"), String::from("{A23C,5.6,C1T,1.2}")),
            ])
        );
    }

    #[test]
    fn parse_attributes_mixed() {
        let mut result = HashMap::new();
        parse_attributes("[&banana;apple=orange]", &mut result);
        assert_eq!(
            result,
            HashMap::from([
                (String::from("banana"), String::from("")),
                (String::from("apple"), String::from("orange")),
            ])
        );
    }

    #[test]
    fn parse_attributes_dups() {
        let mut result = HashMap::new();
        parse_attributes("[&apple=orange;time=2;apple;time=3]", &mut result);
        assert_eq!(
            result,
            HashMap::from([
                (String::from("apple"), String::from("")),
                (String::from("time"), String::from("3")),
            ])
        );
    }
}
