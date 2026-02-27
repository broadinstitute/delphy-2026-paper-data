# Extract date parsing into a library module

## Motivation

`compare_clades.rs` has two date-related functions that are generally
useful for any tree-related tool working with dated tips:

- `parse_tip_date(name: &str) -> Option<f64>`: extracts a date from
  the `name|YYYY-MM-DD` convention and converts to decimal years.
- `decimal_date(year, month, day) -> f64`: converts a calendar date
  to a decimal year, accounting for leap years.

These belong in the library.

## New file: src/dates.rs

Move both functions to a new `dates` module, making them `pub`.
The existing tests for these functions move to the new module as well.

## Changes to src/lib.rs

Add `pub mod dates;`.

## Changes to compare_clades.rs

1. Remove the local `parse_tip_date` and `decimal_date` functions
   and the `// -- Date parsing --` section comment.
2. Add `use tree_ess::dates::{parse_tip_date, decimal_date};`.
3. Remove the 4 date-related tests (`test_parse_tip_date_full`,
   `test_parse_tip_date_partial_month`, `test_parse_tip_date_no_pipe`,
   `test_decimal_date`, `test_decimal_date_leap_year`) since they
   move to the library.

## What is NOT changing

- No behavioral changes.
- The call sites in compare_clades remain the same.
