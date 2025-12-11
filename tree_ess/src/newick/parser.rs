use crate::newick::lexer::Token::Attributes;
use crate::newick::lexer::{
    self, Lexer, LexerError, Token,
    Token::{Colon, Comma, LeftParen, Number, QuotedLabel, RightParen, Semicolon, UnquotedLabel},
};
use crate::newick::parser::ParserError::{Expected, LowerLevel, MultipleRoots, UnbalancedNesting};
use crate::refs::{LightRef, Pool};
use crate::trees::{NodeLike, TreeLike};
use std::collections::HashMap;

#[derive(Debug, PartialEq)]
pub struct NewickTree {
    pub root: LightRef<NewickNode>, // Newick trees always have at least one node
    pub tree_attrs: HashMap<String, String>,
}

impl TreeLike<NewickNode> for NewickTree {
    fn root_ref(&self) -> Option<LightRef<NewickNode>> {
        Some(self.root.clone())
    }
}

#[derive(Debug, PartialEq)]
pub struct NewickNode {
    pub name: String,
    pub node_attrs: HashMap<String, String>,
    pub branch_length: f64,
    pub branch_attrs: HashMap<String, String>,
    pub children: Vec<LightRef<NewickNode>>,
}

impl NodeLike for NewickNode {
    fn child_ref_iter(&self) -> impl DoubleEndedIterator<Item = LightRef<Self>> {
        self.children.iter().map(LightRef::clone)
    }
}

impl NewickTree {
    pub fn new(root: LightRef<NewickNode>) -> Self {
        NewickTree {
            root,
            tree_attrs: HashMap::new(),
        }
    }
    pub fn with_attrs(self, tree_attrs: HashMap<String, String>) -> Self {
        NewickTree { tree_attrs, ..self }
    }
}

impl NewickNode {
    pub fn leaf(name: &str, branch_length: f64) -> NewickNode {
        NewickNode {
            name: String::from(name),
            node_attrs: HashMap::new(),
            branch_length,
            branch_attrs: HashMap::new(),
            children: vec![],
        }
    }

    pub fn inner_node(
        name: &str,
        branch_length: f64,
        children: Vec<LightRef<NewickNode>>,
    ) -> NewickNode {
        NewickNode {
            name: String::from(name),
            node_attrs: HashMap::new(),
            branch_length,
            branch_attrs: HashMap::new(),
            children,
        }
    }

    pub fn with_node_attrs(self, node_attrs: HashMap<String, String>) -> Self {
        NewickNode { node_attrs, ..self }
    }

    pub fn with_branch_attrs(self, branch_attrs: HashMap<String, String>) -> Self {
        NewickNode {
            branch_attrs,
            ..self
        }
    }
}

pub struct Parser<'a> {
    lexer: Lexer<'a>,
    next_token: Result<Token<'a>, LexerError>,
    node_pool: &'a dyn Pool<NewickNode>,
}

#[derive(Debug, PartialEq)]
pub enum ParserError<'a> {
    Expected {
        expected: &'static str,
        got: Token<'a>,
    },
    LowerLevel(LexerError),
    UnbalancedNesting,
    MultipleRoots,
}

impl From<LexerError> for ParserError<'_> {
    fn from(err: LexerError) -> Self {
        LowerLevel(err)
    }
}

impl<'a> Parser<'a> {
    fn new(s: &'a str, pool: &'a dyn Pool<NewickNode>) -> Self {
        let mut lexer = Lexer::new(s);
        let first_token = lexer.next_token();
        Self {
            lexer,
            node_pool: pool,
            next_token: first_token,
        }
    }

    pub fn parse(
        s: &'a str,
        pool: &'a dyn Pool<NewickNode>,
    ) -> Result<NewickTree, ParserError<'a>> {
        let mut parser = Parser::new(s, pool);
        parser.do_parse()
    }

    fn do_parse(&mut self) -> Result<NewickTree, ParserError<'a>> {
        let mut children_stack: Vec<Vec<LightRef<NewickNode>>> = vec![
            vec![], // "Above-root"'s children = roots
            vec![], // partial children of current root
        ];

        let mut tree_attrs: HashMap<String, String> = HashMap::new();
        while let Attributes(tree_attrs_str) = self.peek()? {
            self.advance()?;
            lexer::parse_attributes(&tree_attrs_str, &mut tree_attrs);
        }

        loop {
            // At every iteration of this loop, we're about to start to
            // parse a node that's a child of pending_nodes.top()
            while self.peek()? == LeftParen {
                self.advance()?;
                children_stack.push(Vec::new());
            }

            // We're about to parse a node whose children are in children_stack.last()
            loop {
                let children = children_stack.pop().ok_or(UnbalancedNesting)?;

                // Parse node
                let node = self.parse_rest_of_node(children)?;

                // Add it as a child of parent
                children_stack
                    .last_mut()
                    .ok_or(UnbalancedNesting)?
                    .push(node);

                // Next token determines whether we keep finishing nodes (or the whole tree),
                // or we start a new node
                match self.advance()? {
                    Comma => {
                        // Start a new (currently childless) node that's a child of children_stack.last()
                        children_stack.push(Vec::new());
                        break;
                    }
                    RightParen => {
                        // We're now about to finish parsing a node whose children
                        // are in children_stack.last()
                    }
                    Semicolon => {
                        // The end of the tree!  There should be exactly one item
                        // in `children_stack`, and that is the list of roots of
                        // the current forest.  Though the Newick format could
                        // easily encode forests, we restrict ourselves to single trees
                        let mut roots = children_stack.pop().ok_or(UnbalancedNesting)?;
                        if !children_stack.is_empty() {
                            return Err(UnbalancedNesting);
                        }
                        let root = roots.pop().expect("There's always at least one root!");
                        if !roots.is_empty() {
                            return Err(MultipleRoots);
                        }

                        let mut tree = NewickTree::new(root);
                        tree.tree_attrs = tree_attrs;
                        return Ok(tree);
                    }
                    t => {
                        return Err(Expected {
                            expected: "',', ')' or ';'",
                            got: t,
                        });
                    }
                }
            }
        }
    }

    fn advance(&mut self) -> Result<Token<'a>, LexerError> {
        let cur_token = self.next_token.clone();
        self.next_token = self.lexer.next_token();
        cur_token
    }

    fn peek(&self) -> Result<Token<'a>, LexerError> {
        self.next_token
    }

    fn parse_rest_of_node(
        &mut self,
        children: Vec<LightRef<NewickNode>>,
    ) -> Result<LightRef<NewickNode>, ParserError<'a>> {
        // name=( UNQUOTED_LABEL | QUOTED_LABEL | NUMBER )?
        let name = match self.peek()? {
            UnquotedLabel(label) => {
                self.advance()?;
                String::from(label)
            }
            Number(_, str) => {
                // Use raw text of number in case to avoid information loss, e.g., '00012'
                self.advance()?;
                String::from(str)
            }
            QuotedLabel(label) => {
                self.advance()?;
                lexer::unquote(label).expect("lexer recognized QuotedLabel but unquote failed?")
            }
            _ => String::new(), // No name
        };

        // nodeAttrs=ATTRIBUTES*
        let mut node_attrs = HashMap::new();
        while let Attributes(node_attr_str) = self.peek()? {
            self.advance()?;
            lexer::parse_attributes(node_attr_str, &mut node_attrs);
        }

        // ( ':' branchAttrs=ATTRIBUTES* length=NUMBER? )?
        let mut branch_attrs = HashMap::new();
        let mut branch_length = 0.0;
        if self.peek()? == Colon {
            self.advance()?;
            while let Attributes(branch_attr_str) = self.peek()? {
                self.advance()?;
                lexer::parse_attributes(branch_attr_str, &mut branch_attrs);
            }
            if let Number(num, _) = self.peek()? {
                self.advance()?;
                branch_length = num;
            }
        }

        Ok(self.node_pool.alloc(NewickNode {
            name,
            node_attrs,
            branch_length,
            branch_attrs,
            children,
        }))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refs::TestPool;

    #[test]
    fn empty() {
        let pool = TestPool::new();
        let result = Parser::parse(";", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::leaf("", 0.0))))
        );
    }

    #[test]
    fn single_node_with_unquoted_label() {
        let pool = TestPool::new();
        let result = Parser::parse("A;", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::leaf("A", 0.0))))
        )
    }

    #[test]
    fn single_node_with_quoted_label() {
        let pool = TestPool::new();
        let result = Parser::parse("'A''s node';", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(
                pool.alloc(NewickNode::leaf("A's node", 0.0))
            ))
        )
    }

    #[test]
    fn single_node_with_just_branch_length() {
        let pool = TestPool::new();
        let result = Parser::parse(":123.45;", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::leaf("", 123.45))))
        )
    }

    #[test]
    fn single_node_with_giant_branch_length() {
        let pool = TestPool::new();
        let result = Parser::parse(":1e1000;", &pool); // overflow!

        assert_eq!(
            result,
            Ok(NewickTree::new(
                pool.alloc(NewickNode::leaf("", f64::INFINITY))
            ))
        )
    }

    #[test]
    fn empty_with_empty_child() {
        let pool = TestPool::new();
        let result = Parser::parse("();", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::inner_node(
                "",
                0.0,
                vec![pool.alloc(NewickNode::leaf("", 0.0))]
            ))))
        );
    }

    #[test]
    fn complicated() {
        let pool = TestPool::new();
        let result = Parser::parse("((A:0.1,B:0.2)C:0.3,(D:0.4,E:0.5)F:0.6)R:0.7;", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::inner_node(
                "R",
                0.7,
                vec![
                    pool.alloc(NewickNode::inner_node(
                        "C",
                        0.3,
                        vec![
                            pool.alloc(NewickNode::leaf("A", 0.1)),
                            pool.alloc(NewickNode::leaf("B", 0.2)),
                        ]
                    )),
                    pool.alloc(NewickNode::inner_node(
                        "F",
                        0.6,
                        vec![
                            pool.alloc(NewickNode::leaf("D", 0.4)),
                            pool.alloc(NewickNode::leaf("E", 0.5)),
                        ]
                    )),
                ]
            ))))
        )
    }

    #[test]
    fn with_tree_attrs() {
        let pool = TestPool::new();
        let result = Parser::parse("[&R] [&t=1.2] (A,B);", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::inner_node(
                "",
                0.0,
                vec![
                    pool.alloc(NewickNode::leaf("A", 0.0)),
                    pool.alloc(NewickNode::leaf("B", 0.0)),
                ]
            )))
            .with_attrs(HashMap::from([
                (String::from("R"), String::from("")),
                (String::from("t"), String::from("1.2")),
            ])))
        );
    }

    #[test]
    fn with_overlapping_tree_attrs() {
        let pool = TestPool::new();
        let result = Parser::parse("[&R;t=1.2] [&t=2.4] (A,B);", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::inner_node(
                "",
                0.0,
                vec![
                    pool.alloc(NewickNode::leaf("A", 0.0)),
                    pool.alloc(NewickNode::leaf("B", 0.0)),
                ]
            )))
            .with_attrs(HashMap::from([
                (String::from("R"), String::from("")),
                (String::from("t"), String::from("2.4")),
            ])))
        );
    }

    #[test]
    fn with_overlapping_node_attrs() {
        let pool = TestPool::new();
        let result = Parser::parse("A[&hello;apple] [&hello=world];", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(
                NewickNode::leaf("A", 0.0).with_node_attrs(HashMap::from([
                    (String::from("hello"), String::from("world")),
                    (String::from("apple"), String::from("")),
                ]))
            )))
        );
    }

    #[test]
    fn with_node_attrs() {
        let pool = TestPool::new();
        let result = Parser::parse("(A[&apple=orange],B[&banana]);", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::inner_node(
                "",
                0.0,
                vec![
                    pool.alloc(NewickNode::leaf("A", 0.0).with_node_attrs(HashMap::from([(
                        String::from("apple"),
                        String::from("orange")
                    )]))),
                    pool.alloc(NewickNode::leaf("B", 0.0).with_node_attrs(HashMap::from([(
                        String::from("banana"),
                        String::from("")
                    ),]))),
                ]
            ))))
        );
    }

    #[test]
    fn with_node_and_branch_attrs() {
        let pool = TestPool::new();
        let result = Parser::parse(
            "(A[&apple=orange]:[&muts=None]1.2, B[&banana]:[&muts=Many]);",
            &pool,
        );

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(NewickNode::inner_node(
                "",
                0.0,
                vec![
                    pool.alloc(NewickNode::leaf("A", 1.2).with_node_attrs(HashMap::from([(
                        String::from("apple"),
                        String::from("orange")
                    )])).with_branch_attrs(HashMap::from([
                        (String::from("muts"), String::from("None")),
                    ]))),
                    pool.alloc(
                        NewickNode::leaf("B", 0.0).with_node_attrs(HashMap::from([(
                            String::from("banana"),
                            String::from("")
                        ),])).with_branch_attrs(HashMap::from([
                            (String::from("muts"), String::from("Many")),
                        ]))
                    ),
                ]
            ))))
        );
    }

    #[test]
    fn with_overlapping_branch_attrs() {
        let pool = TestPool::new();
        let result = Parser::parse("A:[&hello;apple] [&hello=world];", &pool);

        assert_eq!(
            result,
            Ok(NewickTree::new(pool.alloc(
                NewickNode::leaf("A", 0.0).with_branch_attrs(HashMap::from([
                    (String::from("hello"), String::from("world")),
                    (String::from("apple"), String::from("")),
                ]))
            )))
        );
    }
}
