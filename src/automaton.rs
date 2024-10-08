use std::cmp;
use std::fmt;

use super::divsufsort::*;
use super::searcher::Searcher;
use super::structs::{ProtoSD, ProtoSDsFamily, RunSettings};
use rayon::prelude::*;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Clone)]
pub struct Segment {
    pub tag: usize,
    pub start: usize,
    pub end: usize,
}

impl fmt::Debug for Segment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "S{{{} -> {}}} ({})",
            self.start,
            self.end,
            self.end - self.start
        )
    }
}
impl Segment {
    fn len(&self) -> usize {
        self.end - self.start
    }
}

#[derive(Debug)]
struct Arm {
    left: Segment,
    right: Segment,
    active: bool,
    dirty: bool,
    gap: usize,
}

enum Operation {
    ExtendArm {
        i: usize,
        l_end: usize,
        r_end: usize,
    },
    NewArm {
        i: usize,
        m_start: usize,
        m_end: usize,
    },
}

#[allow(clippy::too_many_arguments)]
pub fn search_duplications(
    needle: &[u8],
    needle_offset: usize,
    strand: &[u8],
    sa: &[SAIdx],
    searcher: &Searcher,
    progress: &AtomicUsize,
    settings: RunSettings,
) -> Vec<ProtoSDsFamily> {
    fn try_extend_arms(arms: &[Arm], m: &Segment, e: i64, i: usize, ps: usize) -> Operation {
        for (j, a) in arms.iter().enumerate() {
            if a.active
                && d_ss(&a.right, m) < cmp::max(e, (0.1 * a.left.len() as f64) as i64)
                && m.end > a.right.end
            {
                return Operation::ExtendArm {
                    i: j,
                    l_end: i + ps,
                    r_end: m.end,
                };
            }
        }

        Operation::NewArm {
            i,
            m_start: m.start,
            m_end: m.end,
        }
    }

    let mut arms: Vec<Arm> = Vec::new();
    let mut i = 0;
    let mut r = Vec::new();
    let step_size = settings.probe_size / 2;

    if needle.len() < settings.min_duplication_length {
        return Vec::new();
    }

    while i < needle.len() - settings.probe_size - step_size {
        i += step_size;
        progress.store(i, Ordering::Relaxed);

        if needle[i] == b'N' {
            continue;
        }
        let matches: Vec<Segment> = searcher
            .search(strand, sa, &needle[i..i + settings.probe_size])
            .into_iter()
            .filter(|m| m.start != i)
            .filter(|m| {
                if !settings.reverse {
                    m.start > i + needle_offset
                } else {
                    m.start >= needle_offset + needle.len() - i
                }
            })
            .collect();
        if matches.len() > settings.max_cardinality {
            continue;
        }

        // Reset dirty bits of arms
        arms.iter_mut().for_each(|arm| arm.dirty = false);

        let todo = matches
            .par_iter()
            .with_min_len(8)
            .map(|m| {
                try_extend_arms(
                    &arms,
                    m,
                    i64::from(settings.max_gap_size),
                    i,
                    settings.probe_size,
                )
            })
            .collect::<Vec<_>>();

        todo.iter().for_each(|op| {
            if let Operation::ExtendArm { i, l_end, r_end } = op {
                arms[*i].left.end = *l_end;
                arms[*i].right.end = *r_end;
                arms[*i].dirty = true;
                arms[*i].gap = 0;
            }
        });

        todo.iter().for_each(|op| {
            if let Operation::NewArm { i, m_start, m_end } = op {
                arms.push(Arm {
                    left: Segment {
                        start: *i,
                        end: *i + settings.probe_size,
                        tag: 0,
                    },
                    right: Segment {
                        start: *m_start,
                        end: *m_end,
                        tag: 0,
                    },
                    active: true,
                    dirty: false,
                    gap: 0,
                })
            }
        });

        // Update the gaps of non-dirty arms
        arms.iter_mut().filter(|a| !a.dirty).for_each(|a| {
            a.gap += step_size;
            if a.gap as u32 >= settings.max_gap_size {
                a.active = false
            }
        });

        if arms.len() > 200 {
            arms.retain(|a| {
                a.active
                    || a.left.len() >= settings.min_duplication_length
                    || a.right.len() >= settings.min_duplication_length
            });
        }

        // Check if there are still extending arms
        if !arms.is_empty() && arms.iter().all(|a| !a.active) {
            let family: ProtoSDsFamily = arms
                .iter()
                .filter(|a| a.right.len() >= settings.min_duplication_length)
                .map(|a| ProtoSD {
                    left: a.left.start,
                    right: a.right.start,
                    left_length: a.left.len(),
                    right_length: a.right.len(),
                    identity: 0.,
                    reversed: false,
                    complemented: false,
                })
                .collect();
            if !family.is_empty() {
                r.push(family);
            }
            arms.clear();
        }
    }

    r
}

/// Returns the minimal distance between two segments
fn d_ss(a: &Segment, m: &Segment) -> i64 {
    if (m.start >= a.start && m.start <= a.end) || (m.end >= a.start && m.end <= a.end) {
        0
    } else {
        cmp::min(
            (a.start as i64 - m.end as i64).abs(),
            (a.end as i64 - m.start as i64).abs(),
        )
    }
}
