use std::collections::HashMap;

use crate::newick::NewickTree;
use crate::trees::{NodeLike, TraversalAction, TreeLike};

/// Parse a tip date from the `name|YYYY-MM-DD` convention, returning
/// the date as a decimal year.  Returns `None` if the name does not
/// end with a valid full date.
pub fn parse_tip_date(name: &str) -> Option<f64> {
    let date_str = name.rsplit('|').next()?;
    let parts: Vec<&str> = date_str.split('-').collect();
    if parts.len() != 3 {
        return None; // Only full YYYY-MM-DD dates
    }
    let year: i32 = parts[0].parse().ok()?;
    let month: u32 = parts[1].parse().ok()?;
    let day: u32 = parts[2].parse().ok()?;
    if !(1..=12).contains(&month) || !(1..=31).contains(&day) {
        return None;
    }
    Some(decimal_date(year, month, day))
}

/// Convert a calendar date to a decimal year, accounting for leap years.
pub fn decimal_date(year: i32, month: u32, day: u32) -> f64 {
    let is_leap = (year % 4 == 0) && (year % 100 != 0 || year % 400 == 0);
    let days_in_year: f64 = if is_leap { 366.0 } else { 365.0 };
    let days_before_month: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let mut doy = days_before_month[(month - 1) as usize] + day;
    if is_leap && month > 2 {
        doy += 1;
    }
    year as f64 + (doy as f64 - 1.0) / days_in_year
}

/// Use the fixed date of some tip to back out the date of the root.
/// Returns `None` if no tip has an exact date.
///
/// We put up with the possibility of some round-off error in adding and subtracting
/// branch lengths as we traverse the tree.  The working assumption is that most tips
/// have exact dates, so it's very likely that just the initial descent from the root
/// to the first tip suffices to nail down the root date.  If that fails, then likely
/// one of the next few tips visited will have an exact date, thus producing a root date.
/// Hence, the amount of roundoff that we expect in practice is very small.
pub fn find_root_date(
    tree: &NewickTree,
    exact_tip_dates: &HashMap<String, f64>,
) -> Option<f64> {
    let mut root_to_node_dist = 0.0;

    for (action, node_ref) in tree.traversal_iter() {
        let node = node_ref.borrow();
        match action {
            TraversalAction::Enter => {
                if node_ref != tree.root {
                    // root-to-root distance is 0 by definition; ignore branch length then
                    root_to_node_dist += node.branch_length
                };

                if node.is_tip()
                    && let Some(exact_tip_date) = exact_tip_dates.get(node.name.as_str())
                {
                    // Found a tip with an exact date: done!
                    let root_date = exact_tip_date - root_to_node_dist;
                    return Some(root_date);
                }
            }
            TraversalAction::Exit => {
                // If we get here when `node` is the root, then we're going to fail
                // anyway, so there's nothing to gain from subtracting `node.branch_length`
                // only for non-root nodes.
                root_to_node_dist -= node.branch_length;
            }
        }
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_tip_date_full() {
        assert!((parse_tip_date("EBOV-20547|2015-09-07").unwrap() - 2015.6822).abs() < 0.01);
    }

    #[test]
    fn test_parse_tip_date_partial_month() {
        assert_eq!(parse_tip_date("OP415257|9000312|UK|2022-08"), None);
    }

    #[test]
    fn test_parse_tip_date_no_pipe() {
        assert_eq!(parse_tip_date("just_a_name"), None);
    }

    #[test]
    fn test_decimal_date() {
        // Jan 1 of a non-leap year
        assert!((decimal_date(2023, 1, 1) - 2023.0).abs() < 1e-6);
        // Dec 31 of a non-leap year
        assert!((decimal_date(2023, 12, 31) - (2023.0 + 364.0 / 365.0)).abs() < 1e-6);
    }

    #[test]
    fn test_decimal_date_leap_year() {
        // March 1 of a leap year
        assert!((decimal_date(2024, 3, 1) - (2024.0 + 60.0 / 366.0)).abs() < 1e-6);
    }
}
