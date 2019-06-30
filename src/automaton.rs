extern crate rayon;
extern crate indicatif;

use std::cmp;
use std::fmt;

use super::divsufsort64::*;
use super::structs::{RunSettings, ProtoSD, ProtoSDsFamily};
use super::searcher::Searcher;
use std::sync::atomic::{AtomicUsize, Ordering};
use rayon::prelude::*;


#[derive(Clone)]
pub struct Segment {
    pub tag: usize,
    pub start: usize,
    pub end: usize,
}

impl fmt::Debug for Segment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f,
               "S{{{}, {}}} ({})",
               self.start,
               self.end,
               self.end - self.start)
    }
}
impl Segment {
    fn len(&self) -> usize {
        self.end - self.start
    }
}

#[derive(Debug,Clone)]
pub struct ProtoProtoSD {
    bottom: usize,
    top: usize,
    matches: Vec<Segment>,
}



#[derive(Debug)]
struct Arm {
    left: Segment,
    right: Segment,
    family_id: usize,
    active: bool,
    dirty: bool,
    gap: usize,
}

enum Operation {
    ExtendArm {i: usize, l_end: usize, r_end: usize},
    NewArm {i: usize, m_start: usize, m_end: usize},
}

pub fn search_duplications(strand1: &[u8],
                           strand2: &[u8],
                           sa: &[idx],

                           searcher: &Searcher,
                           progress: &AtomicUsize,

                           settings: RunSettings
) -> Vec<ProtoSDsFamily> {
    fn try_extend_arms(arms: &[Arm], m: &Segment, e: i64, i: usize, ps: usize) -> Operation {
        for (j, a) in arms.iter().enumerate() {
            if a.active && d_ss(&a.right, m) < cmp::max(e, (0.1*a.left.len() as f64) as i64) as i64 {
                return Operation::ExtendArm {i: j, l_end: i + ps, r_end: m.end}
            }
        }

        Operation::NewArm {i, m_start: m.start, m_end: m.end}
    }

    let mut arms : Vec<Arm> = Vec::new();
    let mut i = settings.start;
    let mut r = Vec::new();
    let mut current_family_id = 1;
    let step_size = settings.probe_size/2;


    while i < settings.end - settings.probe_size {
        i += step_size; // TODO
        progress.store(cmp::min(i - settings.start, settings.end - settings.start), Ordering::Relaxed);

        if strand1[i] == b'N' { continue }
        let matches: Vec<Segment> = searcher.search(strand2, sa, &strand1[i..i + settings.probe_size])
            .into_iter()
            .filter(|m| m.start != i)
            .filter(|m| if !settings.reverse { m.start > i } else { m.start <= strand2.len() - i })
            .collect();
        if matches.len() > settings.max_cardinality {continue}

        // Reset dirty bits of arms
        arms.iter_mut().for_each(|arm| arm.dirty = false);

        // for m in matches {
        //     // Try to extend existing arms...
        //     if !old_try_extend_arms(&mut arms, &m, settings.max_gap_size as i64, i, settings.probe_size) {
        //         // ...or create new arm
        //         arms.push(Arm{
        //             left: Segment{start: i, end: i + settings.probe_size, tag: 0},
        //             right: Segment{start: m.start, end: m.end, tag: 0},
        //             family_id: current_family_id,
        //             active: true, dirty: false,
        //             gap: 0
        //         })
        //     }
        // }

        let todo = matches
            .par_iter()
            .with_min_len(8)
            .map(|m| try_extend_arms(&arms, m, i64::from(settings.max_gap_size), i, settings.probe_size) )
            .collect::<Vec<_>>();

        todo.iter()
            .for_each(|op| {
                if let Operation::ExtendArm {i, l_end, r_end} = op {
                    arms[*i].left.end = *l_end;
                    arms[*i].right.end = *r_end;
                    arms[*i].dirty = true;
                    arms[*i].gap = 0;
                }
            });

        todo.iter()
            .for_each(|op| {
                if let Operation::NewArm {i, m_start, m_end} = op {
                    arms.push(Arm{
                        left: Segment{start: *i, end: *i + settings.probe_size, tag: 0},
                        right: Segment{start: *m_start, end: *m_end, tag: 0},
                        family_id: current_family_id,
                        active: true, dirty: false,
                        gap: 0
                    })
                }
            });

        // Update the gaps of non-dirty arms
        arms.iter_mut()
            .filter(|a| !a.dirty)
            .for_each(|a|{
                a.gap += step_size;
                if a.gap as u32 >= settings.max_gap_size { a.active = false }
            });

        if arms.len() > 200 {
            arms.retain(|a| {
                a.active || a.left.len() >= settings.min_duplication_length // TODO
            });
        }

        // Check if there are still extending arms
        if !arms.is_empty() && arms.iter().all(|a| !a.active) {
            let family: ProtoSDsFamily = arms.iter()
                .filter(|a| a.right.len() >= settings.min_duplication_length)
                .map(|a| {
                    ProtoSD {
                        left: a.left.start,
                        right: a.right.start,
                        length: a.left.len(),
                        identity: 0.,
                        reversed: false,
                        complemented: false,
                    }}) .collect();
            if !family.is_empty() { r.push(family); }
            arms.clear();

            current_family_id += 1;
        }
    }

    r
}



fn d_ss(a: &Segment, m: &Segment) -> i64 {
    if (m.start >= a.start && m.start <= a.end)
        || (m.end >= a.start && m.end <= a.end) {
            0
        } else {
            cmp::min((a.start as i64 - m.end as i64).abs(),
                     (a.end as i64 - m.start as i64).abs())
        }
}
