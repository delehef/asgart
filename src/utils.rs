use std::cmp;
use std::fmt;
use std::mem;
use std::collections::HashMap;

use super::divsufsort64::*;
use super::structs::SD;

const MIN_DUPLICATION_SIZE: usize = 1000;

pub struct Searcher {
    cache: HashMap<u64, (usize, usize)>,
}
impl Searcher {
    pub fn new(dna: &[u8], sa: &[idx]) -> Searcher {
        let mut s = Searcher {
            cache: HashMap::new()
        };

        let ALPHABET = [b'A', b'T', b'G', b'C', b'N'];
        unsafe{
        for a in ALPHABET.iter() {
            for b in ALPHABET.iter() {
                for c in ALPHABET.iter() {
                    for d in ALPHABET.iter() {
                        for e in ALPHABET.iter() {
                            for f in ALPHABET.iter() {
                                for g in ALPHABET.iter() {
                                    for h in ALPHABET.iter() {
                                        let index = mem::transmute::<[u8; 8], u64>([*a, *b, *c, *d, 
                                                                                   *e, *f, *g, 
                                                                                   *h]);
                                        let p = vec![*a, *b, *c, *d, *e, *f, *g, *h];
                                        let (start, count) : (usize, usize) = {
                                            let mut out = 0;
                                            let count = sa_searchb64(
                                                dna.as_ptr(), dna.len() as i64,
                                                p.as_ptr(), p.len() as i64,
                                                sa.as_ptr(), sa.len() as i64,
                                                &mut out,
                                                0, sa.len() as i64
                                                );
                                            (out as usize, count as usize)
                                        };
                                        s.cache.insert(index, (start, start+count));
                                        // println!("Caching {}", index);
                                        // if count != 0 {
                                        //     println!("Caching {}", 
                                        //              String::from_utf8(p.clone()).unwrap());
                                        //     println!("{}\t{}", start, String::from_utf8(
                                        //             dna[sa[start] as usize..(sa[start]+8) as 
                                        //             usize].to_vec()
                                        //             ).unwrap());
                                        //     println!("...");
                                        //     println!("{}\t{}", start+count-1, String::from_utf8(
                                        //             dna[sa[start+count-1] as 
                                        //             usize..(sa[start+count-1]+8) as 
                                        //             usize].to_vec()
                                        //             ).unwrap());
                                        //     println!("{}\t{}", start+count, String::from_utf8(
                                        //             dna[sa[start+count] as 
                                        //             usize..(sa[start+count]+8) as 
                                        //             usize].to_vec()
                                        //             ).unwrap());

                                        //     let mut pp = p.clone();
                                        //     pp.extend(vec![b'T', b'T', b'G', b'C']);

                                        //     let mut out = 0;
                                        //     let cc = sa_search64(
                                        //         dna.as_ptr(), dna.len() as i64,
                                        //         pp.as_ptr(), pp.len() as i64,
                                        //         sa.as_ptr(), sa.len() as i64,
                                        //         &mut out
                                        //         );
                                        //     println!("First try: {} {}", out, cc);

                                        //     let cc2 = sa_searchb64(
                                        //         dna.as_ptr(), dna.len() as i64,
                                        //         pp.as_ptr(), pp.len() as i64,
                                        //         sa.as_ptr(), sa.len() as i64,
                                        //         &mut out,
                                        //         start as idx, (start+count) as idx
                                        //         );
                                        //     println!("Second try: {} {}", out, cc2);
                                        //     println!("");
                                        //     println!("");
                                        // }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        }

        println!("DONE");
        s
    }

    pub fn search(&self, dna: &[u8], sa: &[idx], pattern: &[u8]) -> Vec<Segment> {
        let index = unsafe {
            mem::transmute::<[u8; 8], u64>(
                [pattern[0], pattern[1], pattern[2], pattern[3], pattern[4], pattern[5], 
                pattern[6], pattern[7]]
                )
        };

        let mut out = 0;
        if !self.cache.contains_key(&index) {
            panic!("Unable to find {} ({})", index, 
                   String::from_utf8(pattern[0..8].to_vec()).unwrap());
        }

        let (lstart, rstart) = self.cache[&index];
        let count = unsafe {
            sa_searchb64(dna.as_ptr(), dna.len() as i64,
            pattern.as_ptr(), pattern.len() as i64,
            sa.as_ptr(), sa.len() as i64, &mut out,
            lstart as idx, rstart as idx)
        };

        let mut rr = Vec::new();
        for i in 0..count {
            let start = sa[(out+i) as usize];
            rr.push(Segment{tag: 0, start: start as usize, end: start as usize + pattern.len()});
        }
        rr
    }
}

#[derive(Clone)]
pub struct Segment {
    pub tag: usize,
    pub start: usize,
    pub end: usize,
}

impl fmt::Debug for Segment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "S{{{}, {}}} ({})", self.start, self.end, self.end - self.start)
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

fn make_duplications(psd: ProtoSD, strand1: &[u8], max_hole_size: u32) -> Vec<SD> {
    let right_segments = make_right_arms(&psd, max_hole_size);
    let mut r = Vec::new();

    for right_segment in right_segments {
        // Ignore too small SD
        let size = right_segment.end - right_segment.start;
        if size < MIN_DUPLICATION_SIZE {continue;}

        // Ignore N-dominant SD
        if *(&strand1[psd.bottom..psd.top].iter().filter(|c| **c == b'n' || **c == b'N').count()) as f32> 0.1*(psd.top - psd.bottom) as f32 {continue;}

        // Ignore overlapping arms
        if ((right_segment.start as i64 - psd.bottom as i64).abs()
            <= (psd.top - psd.bottom) as i64)
            || (right_segment.start == psd.bottom)
            {continue;}

        r.push(SD{
            left: cmp::max(psd.bottom, right_segment.tag),
            right: right_segment.start,
            size: size,
            identity: 0.0
        });
    }

    r
}

pub fn search_duplications(strand1: &[u8], strand2: &[u8], sa: &[idx], start: usize, end: usize,
                           probe_size: usize, max_gap_size: u32,
                           interlaced: bool, searcher: &Searcher) -> Vec<SD> {
    let mut r = Vec::new();

    let mut i = start;
    let mut gap = 0;
    let mut current_start = 0;
    let mut current_segments = Vec::new();
    let mut state = SearchState::Start;

    loop {
        match state {
            SearchState::Start => {
                if i+probe_size >= end-1 {
                    break;
                }
                gap = 0;
                current_start = i;
                if strand1[i] == b'N' || strand1[i] == b'n' {
                    i += 1;
                    state = SearchState::Start;
                } else {
                    current_segments = searcher.search(strand2, sa, &strand1[i..i+probe_size])
                        .into_iter().filter(|x| x.start != i).collect();
                    if current_segments.is_empty() {
                        i += 1;
                        state = SearchState::Start;
                    } else {
                        state = SearchState::Grow;
                    }
                }

            },
            SearchState::Grow => {
                i += probe_size;
                if strand1[i] == b'N' || strand1[i] == b'n' {
                    state = SearchState::SparseGrow;
                } else {
                    let new_matches: Vec<Segment> = searcher.search(strand2, sa, 
                                                                    &strand1[i..cmp::min(i+probe_size, 
                                                                                         strand1.len()-1)], 
                                                                    )
                        .into_iter().filter(|x| x.start != i).collect();

                    if i >= strand1.len() - probe_size {
                        state = SearchState::Proto;
                    } else if segments_to_segments_distance(&new_matches, &current_segments) <= max_gap_size {
                        if interlaced {
                            append_merge_segments(&mut current_segments, &new_matches,
                                                  max_gap_size, i);
                        } else {
                            merge_or_drop_segments(&mut current_segments, &new_matches, max_gap_size);
                        }

                        state = SearchState::Grow;
                    } else {
                        state = SearchState::SparseGrow;
                    }
                }
            },
            SearchState::SparseGrow => {
                i += 1;
                gap += 1;
                if (gap > max_gap_size) || (i >= strand1.len() - probe_size) {
                    state = SearchState::Proto;
                } else if strand1[i] != b'N' && strand1[i] != b'n' {
                    let new_matches = searcher.search(strand2, sa, &strand1[i..i+probe_size]);
                    if segments_to_segments_distance(&new_matches, &current_segments) <= max_gap_size {
                        if interlaced {
                            append_merge_segments(&mut current_segments, &new_matches,
                                                  max_gap_size, i);
                        } else {
                            merge_or_drop_segments(&mut current_segments, &new_matches, max_gap_size);
                        }
                        gap = 0;
                        state = SearchState::Grow;
                    } else {
                        state = SearchState::SparseGrow;
                    }
                }
            },
            SearchState::Proto => {
                if i - current_start >= MIN_DUPLICATION_SIZE {
                    let psd = ProtoSD {
                        bottom: current_start,
                        top: i,
                        matches: current_segments.clone()
                    };
                    let mut result = make_duplications(psd, strand1, max_gap_size);
                    r.append(&mut result);
                }
                state = SearchState::Start;
            },
        }
    }

    r
}

pub fn translate_nucleotide(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' | _ => b'N',
    }
}

pub fn translated(text: &[u8]) -> Vec<u8> {
    text.iter().map(|x| translate_nucleotide(*x)).collect::<Vec<u8>>()
}

fn merge_segments_with_delta(_segments: &[Segment], delta: u64) -> Vec<Segment> {
    let mut segments = _segments.to_vec();
    let mut r = Vec::new();
    segments.sort_by(|x, y| x.start.cmp(&y.start));

    r.push(segments[0].clone());
    for current_segment in segments.iter().skip(1) {
        if current_segment.start - r.last().unwrap().end <= delta as usize{
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
            _originals.push(Segment{tag: i, start: n.start, end: n.end});
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
    let mut distance = 250000000;
    for other in others {
        for segment in segments {
            if other.start >= segment.start && other.start <= segment.end {
                return 0
            } else {
                let local_distance = cmp::min((segment.start as i64 - other.end as i64).abs(), (segment.end as i64 - other.start as i64).abs());
                if local_distance < distance { distance = local_distance }
            }
        }
    }

    distance as u32
}

