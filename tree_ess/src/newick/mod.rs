//! Parsing of Newick-format trees
//!
//! The [Newick format](https://en.wikipedia.org/wiki/Newick_format) is a compact and textual
//! representation for trees.  It lists nodes as they would be visited in a post-order traversal,
//! and uses nested parentheses to encode inner nodes with children.  Each node has an optional
//! label and a distance to its parent node (in externally defined units).
//!
//! A typical Newick tree looks like this:
//! ```text
//! (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)
//! ```
//!
//! The exact Newick tree grammar understood by Delphy, inspired in turn by what is produced by
//! BEAST and its associated tooling (notably TreeAnnotator) is not formally specified anywhere.
//! The closest things to a spec are these three web pages:
//!
//! - <https://phylipweb.github.io/phylip/newicktree.html>
//! - <https://phylipweb.github.io/phylip/newick_doc.html>
//! - <https://beast.community/nexus_metacomments>
//!
//! The grammar below, written in ANTLRv4 format (loosely based
//! on [the example Newick grammar bundled with ANTLRv4](https://github.com/antlr/grammars-v4/blob/master/newick/newick.g4))
//! formally specifies what we recognize:
//!
//! ```text
//! grammar Newick;
//!
//! tree_
//!    : treeAttrs=ATTRIBUTES* root=node ';'
//!    ;
//!
//! node
//!    : ( '(' children+=node (',' children+=node)* ')' )?
//!      name=( UNQUOTED_LABEL | QUOTED_LABEL | NUMBER )?
//!      nodeAttrs=ATTRIBUTES*
//!      ( ':' branchAttrs=ATTRIBUTES* length=NUMBER )?
//!    ;
//!
//! WS             : [ \r\n\t]+ -> skip ;
//! // ATTRIBUTES comes before COMMENT so that "[&hello]" is recognized
//! // as an ATTRIBUTES instead of a COMMENT
//! ATTRIBUTES     : '[&' .*? ']' ;
//! COMMENT        : '[' .*? ']' -> skip ;
//! NUMBER         : [+-]? [0-9]+ ('.' [0-9]+)? ([eE] [+-]? [0-9]+)? ;
//! // Note: unlike spec in https://phylipweb.github.io/phylip/newick_doc.html,
//! // an unquoted label cannot begin with a digit, since that makes distinguishing
//! // NUMBER and UNQUOTED_LABEL tokens impossible.  However, node names are allowed
//! // to be numbers
//! UNQUOTED_LABEL : (~[ \r\n\t()[\]':;,])+ ;
//! QUOTED_LABEL   : '\'' ([^'] | '\'\'')* '\'';
//! ```
//!
//! ## ANTLRv4 Newick grammar license
//!
//! ```text
//! BSD License
//! Copyright (c) 2020, Tom Everett
//! All rights reserved.
//! Redistribution and use in source and binary forms, with or without
//! modification, are permitted provided that the following conditions
//! are met:
//! 1. Redistributions of source code must retain the above copyright
//!    notice, this list of conditions and the following disclaimer.
//! 2. Redistributions in binary form must reproduce the above copyright
//!    notice, this list of conditions and the following disclaimer in the
//!    documentation and/or other materials provided with the distribution.
//! 3. Neither the name of Tom Everett nor the names of its contributors
//!    may be used to endorse or promote products derived from this software
//!    without specific prior written permission.
//! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//! 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
//! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
//! HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
//! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
//! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
//! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//! ```

pub mod lexer;
pub mod parser;

pub use lexer::{Lexer, LexerError};
pub use parser::{Parser, ParserError, NewickNode, NewickTree};