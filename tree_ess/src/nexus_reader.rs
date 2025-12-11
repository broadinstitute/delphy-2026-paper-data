//! Reading of NEXUS files recording posterior trees as produced by BEAST runs.
//!
//! This Nexus reader is *extremely* fragile and not very efficient, but good enough to parse
//! the trees written by Delphy into a `.trees` file and at least one actual BEAST run.

use crate::newick::{NewickNode, NewickTree, Parser};
use crate::refs::Pool;
use crate::trees::TreeLike;
use std::io::{BufRead, BufReader, Read};
use thiserror::Error;

#[derive(Debug, Error)]
pub enum NexusParseError {
    #[error("Unexpected line")]
    UnexpectedLine,

    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    #[error("Invalid number")]
    InvalidNumber,

    #[error("Invalid tree line")]
    InvalidTreeLine,

    #[error("Invalid tip index")]
    InvalidTipIndex,

    #[error("Invalid tree: {0}")]
    InvalidTree(String),
}

struct LineReader<R> {
    reader: R,
    cur_line: String,
}

impl<R: BufRead> LineReader<R> {
    fn new(reader: R) -> Self {
        Self {
            reader,
            cur_line: String::new(),
        }
    }
    fn next_line(&mut self) -> Result<&String, std::io::Error> {
        self.cur_line.clear();
        self.reader.read_line(&mut self.cur_line)?;
        Ok(&self.cur_line)
    }
}

pub struct NexusReader;

impl NexusReader {
    pub fn parse(
        reader: impl Read,
        min_state: u64,
        every: u64,
        pool: &dyn Pool<NewickNode>,
    ) -> Result<Vec<(u64, NewickTree)>, NexusParseError> {
        assert!(every >= 1);

        let reader = BufReader::new(reader);
        let mut line_reader = LineReader::new(reader);

        let mut result = Vec::new();

        fn check(condition: bool) -> Result<(), NexusParseError> {
            if condition {
                Ok(())
            } else {
                Err(NexusParseError::UnexpectedLine)
            }
        }

        check(line_reader.next_line()?.trim() == "#NEXUS")?;
        loop {
            if !line_reader.next_line()?.starts_with('#') {
                break;
            }
        }
        check(line_reader.cur_line.trim() == "")?;
        check(line_reader.next_line()?.trim() == "Begin taxa;")?;
        let ntax_line = line_reader.next_line()?;
        check(ntax_line.trim().starts_with("Dimensions ntax="))?;
        check(ntax_line.trim().ends_with(";"))?;
        let ntax_start = ntax_line.find("=").unwrap() + 1;
        let ntax_end = ntax_line.find(";").unwrap();
        let ntax = ntax_line[ntax_start..ntax_end]
            .parse::<i64>()
            .or(Err(NexusParseError::InvalidNumber))?;
        check(line_reader.next_line()?.trim() == "Taxlabels")?;

        let mut tip_names = vec![];
        tip_names.push(String::from("Tip indices are 1-based, not 0-based!"));
        for _tip_index in 1..=ntax {
            tip_names.push(line_reader.next_line()?.trim().to_string());
        }

        check(line_reader.next_line()?.trim() == ";")?;
        check(line_reader.next_line()?.trim() == "End;")?;

        while line_reader.next_line()?.trim() == "" {}

        check(line_reader.cur_line.trim() == "Begin trees;")?;
        check(line_reader.next_line()?.trim() == "Translate")?;
        // Skip names, we already have them
        for _tip_index in 1..=ntax {
            line_reader.next_line()?;
        }
        check(line_reader.next_line()?.trim() == ";")?;

        loop {
            let line = line_reader.next_line()?;

            if line.trim() == "End;" {
                break;
            }

            check(line.starts_with("tree STATE_"))?;
            let state_start = line.find('_').unwrap() + 1;
            let state_end = line
                .char_indices()
                .skip(state_start) // Works because prefix is ASCII
                .skip_while(|(_, c)| c.is_ascii_digit())
                .next()
                .map(|(i, _)| i)
                .ok_or(NexusParseError::InvalidNumber)?;
            let state = line[state_start..state_end]
                .parse::<u64>()
                .or(Err(NexusParseError::InvalidNumber))?;

            if state < min_state {
                eprintln!("Ignoring state {state} < {min_state} (burn-in)");
            } else if (state % every) != 0 {
                eprintln!("Ignoring state {state} (not a multiple of {every})");
            } else {
                eprintln!("Reading state {state}");

                let newick_start = line.find(" = ").ok_or(NexusParseError::InvalidTreeLine)? + 3;
                let newick_str = line[newick_start..].trim();

                let tree = Parser::parse(newick_str, pool)
                    .map_err(|e| NexusParseError::InvalidTree(format!("{:?}", e)))?;
                for node_ref in tree.iter() {
                    if node_ref.borrow().children.is_empty() {
                        // Translate name index to real name
                        let tip_index = node_ref
                            .borrow()
                            .name
                            .parse::<usize>()
                            .map_err(|_| NexusParseError::InvalidNumber)?;
                        if tip_index >= tip_names.len() {
                            return Err(NexusParseError::InvalidTipIndex);
                        }
                        node_ref.borrow_mut().name = tip_names[tip_index].clone();
                    }
                }
                result.push((state, tree));
            }
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refs::TestPool;
    use indoc::indoc;
    use std::io::Cursor;

    #[test]
    fn simple() {
        let test_contents = indoc! {"\
           #NEXUS
           # A two-line...
           # ...comment

           Begin taxa;
               Dimensions ntax=2;
                   Taxlabels
                       A
                       B
                       ;
           End;
           Begin trees;
               Translate
                   1 A,
                   2 B
           ;
           tree STATE_100 = (1:0.1,2:0.2):0;
           tree STATE_200 = (2:0.3,1:0.4):0;
           End;
        "};

        let pool = TestPool::new();
        let results = NexusReader::parse(Cursor::new(test_contents), 0, 1, &pool);

        assert_eq!(
            results.map_err(|err| err.to_string()),
            Ok(vec![
                (
                    100,
                    NewickTree::new(pool.alloc(NewickNode::inner_node(
                        "",
                        0.0,
                        vec![
                            pool.alloc(NewickNode::leaf("A", 0.1)),
                            pool.alloc(NewickNode::leaf("B", 0.2))
                        ]
                    )))
                ),
                (
                    200,
                    NewickTree::new(pool.alloc(NewickNode::inner_node(
                        "",
                        0.0,
                        vec![
                            pool.alloc(NewickNode::leaf("B", 0.3)),
                            pool.alloc(NewickNode::leaf("A", 0.4))
                        ]
                    )))
                )
            ])
        );
    }
}
