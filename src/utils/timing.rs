use std::cell::RefCell;
use std::collections::BTreeMap;
use std::time::{Duration, Instant};

pub struct Timing {
    pub enabled: bool,
    start: Instant,
    steps: RefCell<BTreeMap<&'static str, Duration>>,
    metrics: RefCell<BTreeMap<&'static str, u128>>,
}

pub struct StepGuard<'a> {
    steps: &'a RefCell<BTreeMap<&'static str, Duration>>,
    label: &'static str,
    began: Instant,
    enabled: bool,
}

impl Timing {
    pub fn new(enabled: bool) -> Self {
        Self {
            enabled,
            start: Instant::now(),
            steps: RefCell::new(BTreeMap::new()),
            metrics: RefCell::new(BTreeMap::new()),
        }
    }
    #[inline]
    pub fn scoped<'a>(&'a self, label: &'static str) -> StepGuard<'a> {
        StepGuard {
            steps: &self.steps,
            label,
            began: Instant::now(),
            enabled: self.enabled,
        }
    }
    #[inline]
    pub fn bump(&self, key: &'static str, delta: u128) {
        if !self.enabled {
            return;
        }
        *self.metrics.borrow_mut().entry(key).or_insert(0) += delta;
    }
    #[inline]
    pub fn set(&self, key: &'static str, v: u128) {
        if !self.enabled {
            return;
        }
        self.metrics.borrow_mut().insert(key, v);
    }
    pub fn finish(&self, title: &str) {
        if !self.enabled {
            return;
        }
        let total = self.start.elapsed();
        eprintln!("\n==== timing report: {title} ====");
        for (k, d) in self.steps.borrow().iter() {
            eprintln!("{k:28} {:>8.3} ms", d.as_secs_f64() * 1e3);
        }
        eprintln!("{:-<44}", "");
        eprintln!("{:<28} {:>8.3} ms", "total", total.as_secs_f64() * 1e3);
        let metrics = self.metrics.borrow();
        if !metrics.is_empty() {
            eprintln!("\n-- metrics --");
            for (k, v) in metrics.iter() {
                eprintln!("{k:28} {v}");
            }
        }
        eprintln!("===============================\n");
    }
}

impl<'a> Drop for StepGuard<'a> {
    fn drop(&mut self) {
        if !self.enabled {
            return;
        }
        let d = self.began.elapsed();
        let mut map = self.steps.borrow_mut();
        *map.entry(self.label).or_insert(Duration::ZERO) += d;
    }
}
