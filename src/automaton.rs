use std::cmp;
use std::fmt;

use super::divsufsort64::*;
use super::structs::SD;
use super::searcher::Searcher;
use std::sync::atomic::{AtomicUsize, Ordering};


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

enum SearchState {
    Start,
    Grow,
    SparseGrow,
    Proto,
}

#[derive(Debug,Clone)]
pub struct ProtoSD {
    bottom: usize,
    top: usize,
    matches: Vec<Segment>,
}

fn make_right_arms(p: &ProtoSD, max_hole_size: u32) -> Vec<Segment> {
    let matches = p.matches.clone();
    merge_segments_with_delta(&matches, max_hole_size as u64)
}


fn make_duplications(psd: ProtoSD,
                     strand1: &[u8],
                     max_hole_size: u32,
                     min_duplication_size: usize)
                     -> Vec<SD> {
    let right_segments = make_right_arms(&psd, max_hole_size);
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

        r.push(SD {
            left: cmp::max(psd.bottom, right_segment.tag),
            right: right_segment.start,
            length: size,
            identity: 0.0,
            reversed: false,
            translated: false,
        });
    }

    r
}

pub fn search_duplications(strand1: &[u8],
                           strand2: &[u8],
                           sa: &[idx],
                           start: usize,
                           end: usize,
                           probe_size: usize,
                           max_gap_size: u32,
                           min_duplication_size: usize,
                           interlaced: bool,
                           searcher: &Searcher,
                           progress: &AtomicUsize)
                           -> Vec<SD> {
    let mut r = Vec::new();

    let mut i = start;
    let mut gap = 0;
    let mut current_start = 0;
    let mut current_segments = Vec::new();
    let mut state = SearchState::Start;

    loop {
        match state {
            SearchState::Start => {
                if i + probe_size >= end - 1 {
                    break;
                }
                gap = 0;
                current_start = i;
                progress.store(cmp::min(i - start, end - start), Ordering::Relaxed);
                state = if strand1[i] == b'N' || strand1[i] == b'n' {
                    i += 1;
                    SearchState::Start
                } else {
                    current_segments = searcher.search(strand2, sa, &strand1[i..i + probe_size])
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
                i += probe_size;
                progress.store(cmp::min(i - start, end - start), Ordering::Relaxed);
                state = if strand1[i] == b'N' || strand1[i] == b'n' {
                    SearchState::SparseGrow
                } else if i + probe_size > strand1.len() - 1 {
                    SearchState::Proto
                } else {
                    let new_matches: Vec<Segment> =
                        searcher.search(strand2, sa, &strand1[i..i + probe_size])
                            .into_iter()
                            .filter(|x| x.start != i)
                            .collect();

                    if segments_to_segments_distance(&new_matches, &current_segments) <=
                       max_gap_size {
                        if interlaced {
                            append_merge_segments(&mut current_segments,
                                                  &new_matches,
                                                  max_gap_size,
                                                  i);
                        } else {
                            merge_or_drop_segments(&mut current_segments,
                                                   &new_matches,
                                                   max_gap_size);
                        }

                        SearchState::Grow
                    } else {
                        SearchState::SparseGrow
                    }
                }
            }
            SearchState::SparseGrow => {
                i += 1;
                gap += 1;
                progress.store(cmp::min(i - start, end - start), Ordering::Relaxed);
                state = if (gap > max_gap_size) || (i >= strand1.len() - probe_size) {
                    SearchState::Proto
                } else if strand1[i] != b'N' && strand1[i] != b'n' {
                    let new_matches = searcher.search(strand2, sa, &strand1[i..i + probe_size]);
                    if segments_to_segments_distance(&new_matches, &current_segments) <=
                       max_gap_size {
                        if interlaced {
                            append_merge_segments(&mut current_segments,
                                                  &new_matches,
                                                  max_gap_size,
                                                  i);
                        } else {
                            merge_or_drop_segments(&mut current_segments,
                                                   &new_matches,
                                                   max_gap_size);
                        }
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
                if i - current_start >= min_duplication_size {
                    let psd = ProtoSD {
                        bottom: current_start,
                        top: i,
                        matches: current_segments.clone(),
                    };
                    let mut result =
                        make_duplications(psd, strand1, max_gap_size, min_duplication_size);
                    r.append(&mut result);
                }
                state = SearchState::Start;
                // println!("{} - {} -> {}", start-start, end-start, i-start);
                progress.store(cmp::min(i - start, end - start), Ordering::Relaxed);
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
    let mut distance = 4000000000;
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
