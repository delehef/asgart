extern crate rayon;
extern crate indicatif;

use std::cmp;
use std::fmt;

use super::divsufsort64::*;
use super::structs::{RunSettings, ProtoSD, ProtoSDsFamily};
use super::searcher::Searcher;
use std::sync::atomic::{AtomicUsize, Ordering};
use rayon::prelude::*;
use self::indicatif::{ProgressBar, ProgressStyle, HumanDuration};


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

enum SearchState {
    Start,
    Grow,
    SparseGrow,
    Proto,
}

#[derive(Debug,Clone)]
pub struct ProtoProtoSD {
    bottom: usize,
    top: usize,
    matches: Vec<Segment>,
}

fn make_right_arms(p: &ProtoProtoSD, max_hole_size: u32) -> Vec<Segment> {
    let matches = p.matches.clone();
    merge_segments_with_delta(&matches, u64::from(max_hole_size))
}


fn make_duplications(psd: &ProtoProtoSD,
                     strand1: &[u8],
                     max_hole_size: u32,
                     min_duplication_size: usize)
                     -> Vec<ProtoSD> {
    let right_segments = make_right_arms(psd, max_hole_size);
    let mut r = Vec::new();

    for right_segment in right_segments {
        // Ignore too small SD
        let size = right_segment.end - right_segment.start;
        if size < min_duplication_size {
            continue;
        }

        // Ignore N-dominant SD
        if strand1[psd.bottom..psd.top]
            .iter()
            .filter(|c| **c == b'n' || **c == b'N')
            .count() as f32 > 0.1 * (psd.top - psd.bottom) as f32 {
                continue;
            }

        // Ignore overlapping arms
        if ((right_segment.start as i64 - psd.bottom as i64).abs() <=
            (psd.top - psd.bottom) as i64) || (right_segment.start == psd.bottom) {
            continue;
        }

        r.push(ProtoSD {
            left: cmp::max(psd.bottom, right_segment.tag),
            right: right_segment.start,
            length: size,
            identity: 0.0,
            reversed: false,
            complemented: false,
        });
    }

    r
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
            if a.active && d_SS(&a.right, m) < cmp::max(e, (0.1*a.left.len() as f64) as i64) as i64 {
                return Operation::ExtendArm {i: j, l_end: i + ps, r_end: m.end}
            }
        }

        return Operation::NewArm {i: i, m_start: m.start, m_end: m.end}
    }

    fn old_try_extend_arms(arms: &mut [Arm], m: &Segment, e: i64, i: usize, ps: usize) -> bool {
        for a in arms {
            if a.active && d_SS(&a.right, m) < e as i64 {
                a.left.end = i + ps;
                a.right.end = m.end;
                a.dirty = true;
                a.gap = 0;
                return true
            }
        }
        false
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
            .map(|m| try_extend_arms(&arms, m, settings.max_gap_size as i64, i, settings.probe_size) )
            .collect::<Vec<_>>();

        todo.iter()
            .for_each(|op| {
                match op {
                    Operation::ExtendArm {i, l_end, r_end} => {
                        arms[*i].left.end = *l_end;
                        arms[*i].right.end = *r_end;
                        arms[*i].dirty = true;
                        arms[*i].gap = 0;
                    }
                    _ => {}
                }
            });

        todo.iter()
            .for_each(|op| {
                    match op {
                        Operation::NewArm {i, m_start, m_end} => {
                            arms.push(Arm{
                                left: Segment{start: *i, end: *i + settings.probe_size, tag: 0},
                                right: Segment{start: *m_start, end: *m_end, tag: 0},
                                family_id: current_family_id,
                                active: true, dirty: false,
                                gap: 0
                            })
                        }
                        _ => {}
                    }
                });

        // Update the gaps of non-dirty arms
        arms.iter_mut()
            .filter(|a| !a.dirty)
            .for_each(|a|{
                a.gap += step_size;
                if a.gap as u32 >= settings.max_gap_size { a.active = false }
            });

        arms.retain(|a| {
            a.active || (!a.active && a.left.len() >= settings.min_duplication_length) // TODO
        });

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
                }})
            .collect();
        if !family.is_empty() { r.push(family); }
        arms.clear();

        current_family_id += 1;
    }
    }

    r
}



pub fn search_duplications_old(strand1: &[u8],
                               strand2: &[u8],
                               sa: &[idx],

                               searcher: &Searcher,
                               progress: &AtomicUsize,

                               settings: RunSettings
) -> Vec<ProtoSD> {
    let mut r = Vec::new();

    let mut prune = 0;
    let mut i = settings.start;
    let mut gap = 0;
    let mut current_start = 0;
    let mut current_segments = Vec::new();
    let mut state = SearchState::Start;

    loop {
        match state {
            SearchState::Start => {
                if i + settings.probe_size >= settings.end - 1 {
                    break;
                }
                gap = 0;
                prune = 0;
                current_start = i;
                progress.store(cmp::min(i - settings.start, settings.end - settings.start), Ordering::Relaxed);
                state = if strand1[i] == b'N' {
                    i += 1;
                    SearchState::Start
                } else {
                    current_segments = searcher.search(strand2, sa, &strand1[i..i + settings.probe_size])
                        .into_iter()
                        .filter(|x| x.start != i)
                        .collect();
                    if current_segments.is_empty() {
                        i += 1;
                        SearchState::Start
                    } else {
                        SearchState::Grow
                    }
                }
            }
            SearchState::Grow => {
                i += settings.probe_size;
                progress.store(cmp::min(i - settings.start, settings.end - settings.start), Ordering::Relaxed);
                state = if strand1[i] == b'N' || strand1[i] == b'n' {
                    SearchState::SparseGrow
                } else if i + settings.probe_size > strand1.len() - 1 {
                    SearchState::Proto
                } else {
                    // prune += 1;
                    // if prune > 2*(settings.min_duplication_length/settings.probe_size) {
                    //     current_segments.retain(|segment|
                    //                             (segment.end - segment.start) > settings.min_duplication_length
                    //     );
                    //     prune = 0;
                    // }

                    let new_matches: Vec<Segment> =
                        searcher.search(strand2, sa, &strand1[i..i + settings.probe_size])
                        .into_iter()
                        .filter(|x| x.start != i)
                        .collect();

                    if new_matches.len() > settings.max_cardinality {
                        if gap > 0 { gap -= 1 };
                        SearchState::SparseGrow
                    } else {
                        if segments_to_segments_distance(&new_matches, &current_segments) <= settings.max_gap_size {
                            // if settings.interlaced {
                            //     append_merge_segments(&mut current_segments,
                            //                           &new_matches,
                            //                           settings.max_gap_size,
                            //                           i);
                            // } else {
                            merge_or_drop_segments(&mut current_segments,
                                                       &new_matches,
                                                       settings.max_gap_size);
                            // }

                            SearchState::Grow
                        } else if current_segments.is_empty() {
                            SearchState::Start
                        } else {
                            SearchState::SparseGrow
                        }
                    }
                }
            }
            SearchState::SparseGrow => {
                i += 1;
                gap += 1;
                progress.store(cmp::min(i - settings.start, settings.end - settings.start), Ordering::Relaxed);
                state = if (gap > settings.max_gap_size) || (i >= strand1.len() - settings.probe_size) {
                    SearchState::Proto
                } else if strand1[i] != b'N' && strand1[i] != b'n' {
                    let new_matches = searcher.search(strand2, sa, &strand1[i..i + settings.probe_size]);
                    if new_matches.len() < settings.max_cardinality
                        && segments_to_segments_distance(&new_matches, &current_segments) <= settings.max_gap_size {
                            // if settings.interlaced {
                            //     append_merge_segments(&mut current_segments,
                            //                           &new_matches,
                            //                           settings.max_gap_size,
                            //                           i);
                            // } else {
                            merge_or_drop_segments(&mut current_segments,
                                                       &new_matches,
                                                       settings.max_gap_size);
                            // }
                            gap = 0;
                            SearchState::Grow
                        } else {
                            SearchState::SparseGrow
                        }
                } else {
                    SearchState::SparseGrow
                }
            }
            SearchState::Proto => {
                if (i - current_start >= settings.min_duplication_length) && !current_segments.is_empty() {
                    let psd = ProtoProtoSD {
                        bottom: current_start,
                        top: i,
                        matches: current_segments.clone(),
                    };
                    let mut result =
                        make_duplications(&psd, strand1, settings.max_gap_size, settings.min_duplication_length);
                    r.append(&mut result);
                }
                state = SearchState::Start;
                progress.store(cmp::min(i - settings.start, settings.end - settings.start), Ordering::Relaxed);
            }
        }
    }

    r
}

fn merge_segments_with_delta(_segments: &[Segment], delta: u64) -> Vec<Segment> {
    let mut segments = _segments.to_vec();
    let mut r = Vec::new();
    segments.sort_by(|x, y| x.start.cmp(&y.start));

    r.push(segments[0].clone());
    for current_segment in segments.iter().skip(1) {
        if current_segment.start as i64 - r.last().unwrap().end as i64 <= delta as i64 {
            r.last_mut().unwrap().end = current_segment.end;
        } else {
            r.push(current_segment.clone());
        }
    }
    r
}

fn append_merge_segments(_originals: &mut Vec<Segment>, _news: &[Segment], delta: u32, i: usize) {
    for n in _news {
        let mut added = false;
        for o in _originals.iter_mut() {
            if n.start >= o.start && n.start <= (o.end + delta as usize) &&
                n.end > (o.end + delta as usize) {
                    o.end = n.end;
                    added = true;
                    continue;
                }
        }
        if !added {
            _originals.push(Segment {
                tag: i,
                start: n.start,
                end: n.end,
            });
        }
    }
}

fn merge_or_drop_segments(_originals: &mut Vec<Segment>, _news: &[Segment], delta: u32) {
    for n in _news {
        for o in _originals.iter_mut() {
            if n.start >= o.start && n.start <= (o.end + delta as usize) &&
                n.end > (o.end + delta as usize) {
                    o.end = n.end;
                    continue;
                }
        }
    }
}

fn segments_to_segments_distance(segments: &[Segment], others: &[Segment]) -> u32 {
    let mut distance = 4_000_000_000;
    for other in others {
        for segment in segments {
            if other.start >= segment.start && other.start <= segment.end {
                return 0;
            } else {
                let local_distance = cmp::min((segment.start as i64 - other.end as i64).abs(),
                                              (segment.end as i64 - other.start as i64).abs());
                if local_distance < distance {
                    distance = local_distance
                }
            }
        }
    }

    distance as u32
}

fn d_SS(a: &Segment, m: &Segment) -> i64 {
    if (m.start >= a.start && m.start <= a.end)
        || (m.end >= a.start && m.end <= a.end) {
            0
        } else {
            cmp::min((a.start as i64 - m.end as i64).abs(),
                     (a.end as i64 - m.start as i64).abs())
        }
}
