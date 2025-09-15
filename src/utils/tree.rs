use serde::{Deserialize, Serialize};

/// Closed interval on [start, end] for point queries.
/// For range queries we use half-open logic.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Interval<T: Ord + Copy> {
    pub start: T,
    pub end: T,
    pub root_fid: u32,
}

#[derive(Debug, Serialize, Deserialize)]
pub struct IntervalTree<T: Ord + Copy> {
    root: Option<Box<Node<T>>>,
}

#[derive(Debug, Serialize, Deserialize)]
struct Node<T: Ord + Copy> {
    center: T,
    intervals: Vec<Interval<T>>,
    left: Option<Box<Node<T>>>,
    right: Option<Box<Node<T>>>,
}

impl<T> IntervalTree<T>
where
    T: Ord + Copy + Serialize + for<'de> Deserialize<'de>,
{
    /// Build a tree from a list of intervals.
    pub fn new(intervals: Vec<Interval<T>>) -> Self {
        let root = Self::build(intervals);
        Self { root }
    }

    fn build(mut intervals: Vec<Interval<T>>) -> Option<Box<Node<T>>> {
        if intervals.is_empty() {
            return None;
        }

        intervals.sort_by_key(|iv| iv.start);
        let mid = intervals.len() / 2;
        let center = intervals[mid].start;

        let mut left = Vec::new();
        let mut right = Vec::new();
        let mut center_ivs = Vec::new();

        for iv in intervals {
            if iv.end < center {
                left.push(iv);
            } else if iv.start > center {
                right.push(iv);
            } else {
                center_ivs.push(iv);
            }
        }

        Some(Box::new(Node {
            center,
            intervals: center_ivs,
            left: Self::build(left),
            right: Self::build(right),
        }))
    }

    /// Point query: returns all intervals covering `point` (closed semantics on [start, end]).
    pub fn query_point(&self, point: T) -> Vec<&Interval<T>> {
        let mut result = Vec::new();
        Self::query_point_rec(&self.root, point, &mut result);
        result
    }

    fn query_point_rec<'a>(
        node: &'a Option<Box<Node<T>>>,
        point: T,
        result: &mut Vec<&'a Interval<T>>,
    ) {
        if let Some(n) = node {
            for iv in &n.intervals {
                if iv.start <= point && point <= iv.end {
                    result.push(iv);
                }
            }
            if point < n.center {
                Self::query_point_rec(&n.left, point, result);
            } else if point > n.center {
                Self::query_point_rec(&n.right, point, result);
            } else {
                // Equal to center: search both sides
                Self::query_point_rec(&n.left, point, result);
                Self::query_point_rec(&n.right, point, result);
            }
        }
    }

    /// Interval query (half-open semantics): returns intervals `iv` where
    /// `iv.start < end && iv.end > start`.
    pub fn query_interval<'a>(&'a self, start: T, end: T, out: &mut Vec<&'a Interval<T>>) {
        Self::query_interval_rec(&self.root, start, end, out);
    }

    fn query_interval_rec<'a>(
        node: &'a Option<Box<Node<T>>>,
        start: T,
        end: T,
        out: &mut Vec<&'a Interval<T>>,
    ) {
        if let Some(n) = node {
            for iv in &n.intervals {
                if iv.start < end && iv.end > start {
                    out.push(iv);
                }
            }
            if start < n.center {
                Self::query_interval_rec(&n.left, start, end, out);
            }
            if end > n.center {
                Self::query_interval_rec(&n.right, start, end, out);
            }
        }
    }
}