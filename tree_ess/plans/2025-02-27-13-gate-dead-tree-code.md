# Gate dead code in trees.rs with #[cfg(test)]

## Context

`trees.rs` contains several concrete tree implementations that are only
used in the module's own tests, generating compiler warnings.

## Changes to src/trees.rs

### Gate `SimpleNode` / `SimpleTree` with `#[cfg(test)]`

`SimpleNode<T>`, `SimpleTree<T>`, and their `TreeLike` / `NodeLike`
impls (lines 41-61) are only used by the tests in this module.  Gate
them with `#[cfg(test)]`.

### Gate `BinaryNode` / `BinaryTree` / `MyIter` with `#[cfg(test)]`

`BinaryNode<T>`, `BinaryTree<T>`, `MyIter<T>`, and their trait impls
(lines 63-125) are not used anywhere, not even in tests.  Gate them
with `#[cfg(test)]` to preserve the assembly-quality documentation
in the comment block above `BinaryNode`.

### Remove commented-out `VecBackedTree` code

The commented-out `VecBackedTree` block (lines 372-466) is abandoned
exploratory code.  Remove it entirely.

## What is NOT changing

- `TreeLike`, `NodeLike`, `TraversalAction`, and the traversal
  iterators are unchanged — they are used by production code.

## Verification

`cargo test` — all tests should pass, and the dead_code warnings
should disappear.
