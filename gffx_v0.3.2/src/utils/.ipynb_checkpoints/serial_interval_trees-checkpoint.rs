use anyhow::Result;
use bincode2::{deserialize, deserialize_from, serialize, serialize_into};
use serde::{Deserialize, Serialize, de::DeserializeOwned};
use serde_json;
use std::{
    fs::{File, write},
    io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write},
    str,
};

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

impl<T: Ord + Copy + Serialize + for<'de> Deserialize<'de>> IntervalTree<T> {
    pub fn new(intervals: Vec<Interval<T>>) -> Self {
        let root = Self::build(intervals);
        Self { root }
    }

    pub fn save_to_file(&self, path: &str) -> std::io::Result<()> {
        let encoded = serialize(self).expect("Serialization failed");
        write(path, encoded)?;
        Ok(())
    }

    pub fn load_from_file(path: &str) -> std::io::Result<Self> {
        let mut file = File::open(path)?;
        let mut buf = Vec::new();
        file.read_to_end(&mut buf)?;
        let tree: Self = deserialize(&buf).expect("Deserialization failed");
        Ok(tree)
    }

    pub fn query_point(&self, point: T) -> Vec<&Interval<T>> {
        let mut result = Vec::new();
        Self::query_point_rec(&self.root, point, &mut result);
        result
    }

    pub fn query_interval(&self, start: T, end: T) -> Vec<&Interval<T>> {
        let mut result = Vec::new();
        Self::query_interval_rec(&self.root, start, end, &mut result);
        result
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
                // == center，search in both directions
                Self::query_point_rec(&n.left, point, result);
                Self::query_point_rec(&n.right, point, result);
            }
        }
    }

    fn query_interval_rec<'a>(
        node: &'a Option<Box<Node<T>>>,
        start: T,
        end: T,
        result: &mut Vec<&'a Interval<T>>,
    ) {
        if let Some(n) = node {
            for iv in &n.intervals {
                if iv.start <= end && iv.end >= start {
                    result.push(iv);
                }
            }
            if start < n.center {
                Self::query_interval_rec(&n.left, start, end, result);
            }
            if end > n.center {
                Self::query_interval_rec(&n.right, start, end, result);
            }
        }
    }
}

pub fn save_multiple_trees<T: Ord + Copy + Serialize>(
    output_path: &str,
    trees: &[IntervalTree<T>],
) -> Result<Vec<u64>> {
    let file = File::create(output_path)?;
    let mut writer = BufWriter::new(file);
    let mut offsets = Vec::with_capacity(trees.len());

    for tree in trees {
        let pos = writer.stream_position()?;
        offsets.push(pos);

        serialize_into(&mut writer, tree).expect("Failed to serialize tree");
    }

    writer.flush()?;
    Ok(offsets)
}

pub fn write_offsets_to_file(offsets: &[u64], path: &str) -> Result<()> {
    let json = serde_json::to_string_pretty(offsets).expect("Serialize offsets failed");
    Ok(write(path, json)?)
}

pub fn load_multiple_trees<T>(rit_path: &str, rix_path: &str) -> Result<Vec<IntervalTree<T>>>
where
    T: Ord + Copy + DeserializeOwned,
{
    let mut reader = BufReader::new(File::open(rit_path)?);
    let offsets: Vec<u64> = {
        let f = File::open(rix_path)?;
        serde_json::from_reader(f)?
    };

    let mut trees = Vec::with_capacity(offsets.len());
    for &off in &offsets {
        reader.seek(SeekFrom::Start(off))?;
        let tree: IntervalTree<T> = deserialize_from(&mut reader)
            .map_err(|e| anyhow::anyhow!("Deserializing failed: {}", e))?;
        trees.push(tree);
    }
    Ok(trees)
}
