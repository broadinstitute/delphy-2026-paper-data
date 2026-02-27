// TODO: Abstract NodeRef<T> so that we can use something punchier than a LightRef<N>
//       to refer to nodes when needed (e.g., a Vec-backed tree).

use crate::refs::LightRef;

pub trait TreeLike<N: NodeLike>: Sized {
    fn root_ref(&self) -> Option<LightRef<N>>;

    // For when the order doesn't matter: the most efficient order
    fn iter(&self) -> impl Iterator<Item=LightRef<N>> {
        self.any_order_iter()
    }

    // For when the order doesn't matter: the most efficient order
    fn any_order_iter(&self) -> impl Iterator<Item=LightRef<N>> {
        self.preorder_iter()
    }

    fn preorder_iter(&self) -> impl Iterator<Item=LightRef<N>> {
        PreorderIterator::new(self)
    }

    fn postorder_iter(&self) -> impl Iterator<Item=LightRef<N>> {
        PostorderIterator::new(self)
    }

    fn traversal_iter(&self) -> impl Iterator<Item=(TraversalAction, LightRef<N>)> {
        TraversalIterator::new(self)
    }
}

pub trait NodeLike: Sized {
    fn child_ref_iter(&self) -> impl DoubleEndedIterator<Item=LightRef<Self>>;

    fn is_tip(&self) -> bool {
        self.child_ref_iter().next().is_none()
    }
}

// Various traversals
// ==================

// Preorder iteration
// ------------------

struct PreorderIterator<N: NodeLike> {
    stack: Vec<LightRef<N>>,
}

impl<N: NodeLike> PreorderIterator<N>
{
    fn new(tree: &impl TreeLike<N>) -> Self {
        PreorderIterator {
            stack: tree.root_ref().iter().map(LightRef::clone).collect(),
        }
    }
}

impl<N: NodeLike> Iterator for PreorderIterator<N> {
    type Item = LightRef<N>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop() {
            None => None,
            Some(node_ref) => {
                for child_ref in node_ref.borrow().child_ref_iter().rev() {
                    self.stack.push(child_ref.clone());
                }

                Some(node_ref)
            }
        }
    }
}


// Post-order iteration
// --------------------

enum PostorderAction {
    Expand,
    Visit,
}

struct PostorderIterator<N: NodeLike> {
    stack: Vec<(PostorderAction, LightRef<N>)>,
}

impl<N: NodeLike> PostorderIterator<N> {
    fn new(tree: &impl TreeLike<N>) -> Self {
        PostorderIterator {
            stack: tree.root_ref().iter().map(|r| (PostorderAction::Expand, r.clone())).collect(),
        }
    }
}

impl<N: NodeLike> Iterator for PostorderIterator<N> {
    type Item = LightRef<N>;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            match self.stack.pop() {
                None => return None,
                Some((PostorderAction::Visit, node_ref)) => return Some(node_ref),
                Some((PostorderAction::Expand, node_ref)) => {
                    self.stack.push((PostorderAction::Visit, node_ref.clone()));
                    for child_ref in node_ref.borrow().child_ref_iter().rev() {
                        self.stack
                            .push((PostorderAction::Expand, child_ref.clone()));
                    }
                }
            }
        }
    }
}

// Explicit traversals
// -------------------
// These iterate through every Entry and Exit of every node during a tree traversal

#[derive(Debug)]
#[derive(PartialEq)]
pub enum TraversalAction {
    Enter,
    Exit,
}
struct TraversalIterator<N: NodeLike> {
    stack: Vec<(TraversalAction, LightRef<N>)>,
}
impl<N: NodeLike> TraversalIterator<N> {
    fn new(tree: &impl TreeLike<N>) -> Self {
        TraversalIterator {
            stack: tree.root_ref().iter().map(|r| (TraversalAction::Enter, r.clone())).collect(),
        }
    }
}
impl<N: NodeLike> Iterator for TraversalIterator<N> {
    type Item = (TraversalAction, LightRef<N>);

    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop() {
            None => None,
            Some((TraversalAction::Enter, node_ref)) => {
                self.stack.push((TraversalAction::Exit, node_ref.clone()));
                for child_ref in node_ref.borrow().child_ref_iter().rev() {
                    self.stack
                        .push((TraversalAction::Enter, child_ref.clone()));
                }
                Some((TraversalAction::Enter, node_ref))
            }
            Some((TraversalAction::Exit, node_ref)) => {
                Some((TraversalAction::Exit, node_ref))
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::refs::{TestPool, Pool};
    use std::iter;

    // -- N-ary tree (SimpleNode / SimpleTree) --

    #[derive(Debug, PartialEq, Eq)]
    struct SimpleNode<T> {
        value: T,
        child_refs: Vec<LightRef<SimpleNode<T>>>,
    }

    #[derive(Debug, PartialEq, Eq)]
    struct SimpleTree<T> {
        root_ref: Option<LightRef<SimpleNode<T>>>,
    }

    impl<T> TreeLike<SimpleNode<T>> for SimpleTree<T> {
        fn root_ref(&self) -> Option<LightRef<SimpleNode<T>>> {
            self.root_ref.clone()
        }
    }
    impl<T> NodeLike for SimpleNode<T> {
        fn child_ref_iter(&self) -> impl DoubleEndedIterator<Item=LightRef<SimpleNode<T>>> {
            self.child_refs.iter().map(LightRef::clone)
        }
    }

    // -- Binary tree (BinaryNode / BinaryTree) --
    //
    // AMAZINGLY, deduced from staring at long assembly listing, when using LightRef<T> =~ NonNull<T> :
    //
    // (1) child_refs (: Option<(LightRef<BinaryNode<T>>, LightRef<BinaryNode<T>>)>)
    //     really is stored as two pointers, with a NULL on the first pointer signalling
    //     `None`
    // (2) All the child_ref_iter() boilerplate (even child_ref_iter().rev()) gets *completely*
    //     eliminated, and is replaced by testing for a leaf (left == NULL), and if not,
    //     a direct reading of the left and right pointers!
    //
    // Indeed, a preorder iteration through a tree with a simple action on every visit is compiled
    // into nearly optimal assembly of pointer chasing!  (a few needless stores to the stack
    // that don't cause hazards => they go to L1 and are almost hidden by ILP; and a needless test
    // for an empty stack right after pushing left & right pointers, which is almost hidden by
    // what should be perfect branch prediction).
    //
    // Rust + LLVM are extraordinary...

    #[derive(Debug, PartialEq, Eq)]
    struct BinaryNode<T> {
        value: T,
        child_refs: Option<(LightRef<BinaryNode<T>>, LightRef<BinaryNode<T>>)>,
    }

    #[derive(Debug, PartialEq, Eq)]
    struct BinaryTree<T> {
        root_ref: Option<LightRef<BinaryNode<T>>>,
    }

    impl<T> TreeLike<BinaryNode<T>> for BinaryTree<T> {
        fn root_ref(&self) -> Option<LightRef<BinaryNode<T>>> {
            self.root_ref.clone()
        }
    }
    enum MyIter<T> {
        Zero(iter::Empty<T>),
        Two(std::array::IntoIter<T, 2>),
    }
    impl<T> Iterator for MyIter<T> {
        type Item = T;
        fn next(&mut self) -> Option<Self::Item> {
            match self {
                MyIter::Zero(iter) => iter.next(),
                MyIter::Two(iter) => iter.next(),
            }
        }
    }
    impl<T> DoubleEndedIterator for MyIter<T> {
        fn next_back(&mut self) -> Option<Self::Item> {
            match self {
                MyIter::Zero(iter) => iter.next_back(),
                MyIter::Two(iter) => iter.next_back(),
            }
        }
    }
    impl<T> NodeLike for BinaryNode<T> {
        fn child_ref_iter(&self) -> impl DoubleEndedIterator<Item=LightRef<BinaryNode<T>>> {
            match &self.child_refs {
                None => MyIter::Zero(iter::empty()),
                Some((l, r)) => MyIter::Two([l.clone(), r.clone()].into_iter()),
            }
        }
    }

    // -- Tests --

    #[test]
    fn any_order_iter() {
        let pool = TestPool::new();
        let child1_ref = pool.alloc(SimpleNode {
            value: 2,
            child_refs: vec![],
        });
        let child2_ref = pool.alloc(SimpleNode {
            value: 3,
            child_refs: vec![],
        });
        let root_ref = pool.alloc(SimpleNode {
            value: 1,
            child_refs: vec![child1_ref.clone(), child2_ref.clone()],
        });
        let tree = SimpleTree {
            root_ref: Some(root_ref.clone()),
        };

        let mut sum_values = 0;
        for node_ref in tree.any_order_iter() {
            sum_values += node_ref.borrow().value;
        }
        assert_eq!(sum_values, 6);

        drop(tree);

        // Check that these objects can be explicitly freed in this order
        pool.free(root_ref);
        pool.free(child1_ref);
        pool.free(child2_ref);
    }

    #[test]
    fn any_order_iter_mut() {
        let pool = TestPool::new();
        let child1_ref = pool.alloc(SimpleNode {
            value: 2,
            child_refs: vec![],
        });
        let child2_ref = pool.alloc(SimpleNode {
            value: 3,
            child_refs: vec![],
        });
        let root_ref = pool.alloc(SimpleNode {
            value: 1,
            child_refs: vec![child1_ref.clone(), child2_ref.clone()],
        });
        let tree = SimpleTree {
            root_ref: Some(root_ref.clone()),
        };

        for node_ref in tree.any_order_iter() {
            node_ref.borrow_mut().value *= 2;
        }
        let values = tree
            .preorder_iter()
            .map(|node_ref| node_ref.borrow().value)
            .collect::<Vec<i32>>();
        assert_eq!(values, vec![2, 4, 6]);
    }

    #[test]
    fn post_order_iter() {
        let pool = TestPool::new();
        let tree = SimpleTree {
            root_ref: Some(pool.alloc(SimpleNode {
                value: 1,
                child_refs: vec![
                    pool.alloc(SimpleNode {
                        value: 2,
                        child_refs: vec![],
                    }),
                    pool.alloc(SimpleNode {
                        value: 3,
                        child_refs: vec![],
                    }),
                ],
            })),
        };

        let values = tree
            .postorder_iter()
            .map(|node_ref| node_ref.borrow().value)
            .collect::<Vec<i32>>();
        assert_eq!(values, vec![2, 3, 1]);
    }

    #[test]
    fn traversal_iter() {
        let pool = TestPool::new();
        let tree = SimpleTree {
            root_ref: Some(pool.alloc(SimpleNode {
                value: 1,
                child_refs: vec![
                    pool.alloc(SimpleNode {
                        value: 2,
                        child_refs: vec![],
                    }),
                    pool.alloc(SimpleNode {
                        value: 3,
                        child_refs: vec![],
                    }),
                ],
            })),
        };

        let actions: Vec<(TraversalAction, i64)> = tree
            .traversal_iter()
            .map(|(action, node_ref)| (action, node_ref.borrow().value))
            .collect();
        assert_eq!(actions, vec![
            (TraversalAction::Enter, 1),
            (TraversalAction::Enter, 2),
            (TraversalAction::Exit, 2),
            (TraversalAction::Enter, 3),
            (TraversalAction::Exit, 3),
            (TraversalAction::Exit, 1),
        ]);
    }

    // Analogous tests for BinaryNode / BinaryTree

    #[test]
    fn binary_any_order_iter() {
        let pool = TestPool::new();
        let child1_ref = pool.alloc(BinaryNode {
            value: 2,
            child_refs: None,
        });
        let child2_ref = pool.alloc(BinaryNode {
            value: 3,
            child_refs: None,
        });
        let tree = BinaryTree {
            root_ref: Some(pool.alloc(BinaryNode {
                value: 1,
                child_refs: Some((child1_ref, child2_ref)),
            })),
        };

        let mut sum_values = 0;
        for node_ref in tree.any_order_iter() {
            sum_values += node_ref.borrow().value;
        }
        assert_eq!(sum_values, 6);
    }

    #[test]
    fn binary_post_order_iter() {
        let pool = TestPool::new();
        let tree = BinaryTree {
            root_ref: Some(pool.alloc(BinaryNode {
                value: 1,
                child_refs: Some((
                    pool.alloc(BinaryNode { value: 2, child_refs: None }),
                    pool.alloc(BinaryNode { value: 3, child_refs: None }),
                )),
            })),
        };

        let values: Vec<i32> = tree
            .postorder_iter()
            .map(|node_ref| node_ref.borrow().value)
            .collect();
        assert_eq!(values, vec![2, 3, 1]);
    }

    #[test]
    fn binary_traversal_iter() {
        let pool = TestPool::new();
        let tree = BinaryTree {
            root_ref: Some(pool.alloc(BinaryNode {
                value: 1,
                child_refs: Some((
                    pool.alloc(BinaryNode { value: 2, child_refs: None }),
                    pool.alloc(BinaryNode { value: 3, child_refs: None }),
                )),
            })),
        };

        let actions: Vec<(TraversalAction, i32)> = tree
            .traversal_iter()
            .map(|(action, node_ref)| (action, node_ref.borrow().value))
            .collect();
        assert_eq!(actions, vec![
            (TraversalAction::Enter, 1),
            (TraversalAction::Enter, 2),
            (TraversalAction::Exit, 2),
            (TraversalAction::Enter, 3),
            (TraversalAction::Exit, 3),
            (TraversalAction::Exit, 1),
        ]);
    }
}
