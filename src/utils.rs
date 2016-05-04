use std::collections::HashSet;
use std::cmp;
use bio::data_structures::suffix_array::SuffixArray;
use std::io::Write;
use std::fmt;
use bio::alignment::pairwise::*;

const MIN_PALINDROME_SIZE: usize = 1000;
const MAX_ALIGNMENT_SIZE: usize = 100000;


macro_rules! log(
    ($($arg:tt)*) => (
        // match writeln!(&mut ::std::io::stderr(), $($arg)* ) {
        //     Ok(_) => {},
        //     Err(x) => {} //panic!("Unable to write to stderr: {}", x),
        // }
    )
);


#[derive(Debug, Clone)]
pub struct Palindrome {
    pub left: usize,
    pub right: usize,
    pub size: usize,
    pub rate: f32,
}

enum SearchState {
    Start,
    Grow,
    SparseGrow,
    Proto,
}

#[derive(Debug,Clone)]
pub struct ProtoPalindrome {
    bottom: usize,
    top: usize,
    matches: Vec<Segment>,
}

#[derive(Debug,Clone)]
pub enum ProcessingPalindrome {
    Done(Palindrome),
    ForFuzzy{pp: ProtoPalindrome, right: usize},
    ForSW{pp: ProtoPalindrome, left_match: Vec<u8>, right_match: Vec<u8>, right_segment: Segment},
    Empty,
}

fn make_right_arms(p: &ProtoPalindrome, max_hole_size: u32) -> Vec<Segment> {
    let matches = p.matches.clone();
    merge_segments_with_delta(&matches, max_hole_size as u64)
}

fn make_right_arm(p: &ProtoPalindrome, max_hole_size: u32) -> Segment {
    let mut matches = p.matches.clone();
    matches = merge_segments_with_delta(&matches, max_hole_size as u64);

    // Sort by size, descending
    matches.sort_by(|a, b| (b.end - b.start).cmp(&(a.end - a.start)));

    return matches.remove(0);
}

fn make_palindrome(pp: ProtoPalindrome, dna: &[u8], rt_dna: &[u8], max_hole_size: u32) -> Vec<ProcessingPalindrome> {
    let right_segments = make_right_arms(&pp, max_hole_size);
    let mut r = Vec::new();

    for right_segment in right_segments.iter() {
        if right_segment.end - right_segment.start < MIN_PALINDROME_SIZE {continue;}

        // Fetch left and right candidate areas
        let left_match = &dna[pp.bottom..pp.top];
        let right_match = &rt_dna[right_segment.start..right_segment.end];

        r.push(ProcessingPalindrome::Done(Palindrome{left: pp.bottom, right: right_segment.start, size: pp.top - pp.bottom, rate: 0.0}));
        // if right_match.len() > MAX_ALIGNMENT_SIZE || left_match.len() > MAX_ALIGNMENT_SIZE {
        //     r.push(ProcessingPalindrome::ForFuzzy {
        //         pp: pp.clone(),
        //         right: right_segment.start,
        //     });
        // } else {
        //     r.push(ProcessingPalindrome::ForSW {
        //         pp: pp.clone(),
        //         left_match: Vec::from(left_match),
        //         right_match: Vec::from(right_match),
        //         right_segment: right_segment.clone(),
        //     })
        // }
    }

    r
}

pub fn sw_align(p: ProcessingPalindrome) -> ProcessingPalindrome {
    if let ProcessingPalindrome::ForSW{pp, left_match, right_match, right_segment} = p {
        let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
        let mut aligner = Aligner::with_capacity(left_match.len(), right_match.len(), -5, -1, &score);
        let alignment = aligner.local(&left_match, &right_match);
        if alignment.operations.len() < MIN_PALINDROME_SIZE {return ProcessingPalindrome::Empty;};

        let result = Palindrome {
            left: pp.bottom + alignment.xstart,
            right: right_segment.start + alignment.ystart,
            size: alignment.operations.len(),
            rate: 0.0
        };
        println!("Left position:  {}", result.left);
        println!("Right position: {}", result.right);
        println!("Size:           {}", result.size);
        println!("Alignment: \n{}\n", alignment.pretty(&left_match, &right_match));

        ProcessingPalindrome::Done(result)
    } else {
        p
    }
}

pub fn fuzzy_align(dna: &[u8], rt_dna: &[u8], p: ProcessingPalindrome, max_hole_size: u32) -> ProcessingPalindrome {
    if let ProcessingPalindrome::ForFuzzy{pp, ..} = p {
        let moar = expand_nw(dna, rt_dna, pp.bottom, make_right_arm(&pp, max_hole_size).start);
        ProcessingPalindrome::Done(moar)
    } else {
        p
    }
}

pub fn make_palindromes(dna: &[u8], rt_dna: &[u8], sa: &SuffixArray, start: usize, end: usize,
                        candidate_size: usize, max_hole_size: u32) -> Vec<ProcessingPalindrome> {
    let mut r = Vec::new();

    let mut i = start;
    let mut hole = 0;
    let mut current_start = 0;
    let mut current_segments = Vec::new();
    let mut state = SearchState::Start;

    loop {
        match state {
            SearchState::Start => {
                if i+candidate_size >= end-1 {
                    break;
                }
                hole = 0;
                current_start = i;
                current_segments = Vec::new();
                let mut new_segments = search(&rt_dna, &sa, &dna[i..i+candidate_size], candidate_size);
                if new_segments.len() > 0 {
                    current_segments = merge_segments(&mut new_segments);
                    state = SearchState::Grow;
                } else {
                    i += 1;
                    state = SearchState::Start;
                }

            },
            SearchState::Grow => {
                i += candidate_size/2;
                let mut new_matches = search(&rt_dna, &sa, &dna[i..i+candidate_size], candidate_size);

                if i >= dna.len() - candidate_size {
                    state = SearchState::Proto;
                } else if segments_to_segments_distance(&new_matches, &current_segments) <= max_hole_size {
                    current_segments.append(&mut new_matches);
                    current_segments = merge_segments(&mut current_segments);
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::SparseGrow => {
                i += 1;
                hole += 1;
                let mut new_matches = search(&rt_dna, &sa, &dna[i..i+candidate_size], candidate_size);

                if hole > max_hole_size {
                    state = SearchState::Proto;
                } else if i >= dna.len() - candidate_size {
                    state = SearchState::Proto;
                } else if segments_to_segments_distance(&new_matches, &current_segments) <= max_hole_size {
                    current_segments.append(&mut new_matches);
                    current_segments = merge_segments(&mut current_segments);
                    hole = 0;
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::Proto => {
                if i - current_start >= MIN_PALINDROME_SIZE {
                    let pp = ProtoPalindrome {
                        bottom: current_start,
                        top: i,
                        matches: current_segments.clone()
                    };
                    let mut result = make_palindrome(pp, dna, rt_dna, max_hole_size);
                    r.append(&mut result);
                }
                state = SearchState::Start;
            },
        }
    }

    r
}


pub struct AlignmentResult {
    pub errors: u32,
    pub ins_la: u32,
    pub ins_ra: u32,
}


pub fn needleman_wunsch(q1: &[u8], q2: &[u8]) -> AlignmentResult {
    const DELETE_SCORE: i64 = -5;
    const MISMATCH_SCORE: i64 = -3;
    const MATCH_SCORE: i64 = 1;

    assert!(q1.len() == q2.len());
    let size: usize = q1.len();

    let mut F = vec![0 as i64; size * size];

    for i in 0..size {
        F[i] = DELETE_SCORE*i as i64;
    }

    for j in 0..size {
        F[j*size as usize] = DELETE_SCORE*j as i64;
    }

    for i in 1..size {
        for j in 1..size {
            let _match = F[(i-1)+(j-1)*size] + if q1[i] == q2[j] {MATCH_SCORE} else {MISMATCH_SCORE};
            let del = F[(i-1)+(j)*size] + DELETE_SCORE;
            let ins = F[i+(j-1)*size] + DELETE_SCORE;
            F[i+j*size] = cmp::max(_match, cmp::max(del, ins));
        }
    }

    let mut i = size-1 as usize;
    let mut j = size-1 as usize;

    let mut result = AlignmentResult {errors: 0, ins_la: 0, ins_ra: 0};
    while i>1 && j>1 {
        if i>0 && j>0
            && (F[i+j*size] == F[(i-1)+(j-1)*size] + if q1[i] == q2[j] {MATCH_SCORE} else {MISMATCH_SCORE})
            {
                if q1[i] != q2[j] {
                    result.errors += 1;
                }
                i -= 1;
                j -= 1;
            } else if i>0 && F[i+j*size] == F[(i-1)+j*size] + DELETE_SCORE {
                result.ins_la += 1;
                i -= 1;
            } else if j>0 && F[i+j*size] == F[i+(j-1)*size] + DELETE_SCORE {
                result.ins_ra += 1;
                j -= 1;
            }
    }
    result
}


pub fn strcmp(s1: &[u8], s2: &[u8]) -> i32 {
    assert!(s1.len() == s2.len());
    let n = s2.len();

    for i in 0..n {
        if s1[i] != s2[i] {
            return s1[i] as i32 - s2[i] as i32;
        }
    }
    0
}


pub fn translate_nucleotide(n: u8) -> u8 {
    if n == 'A' as u8 { 'T' as u8}
    else if n == 'T' as u8 { 'A' as u8 }
    else if n == 'G' as u8 { 'C' as u8 }
    else if n == 'C' as u8 { 'G' as u8 }
    else if n == 'N' as u8 { 'N' as u8 }
    else { panic!("Not DNA: {}", n as u8); }
}


pub fn translated(text: &[u8]) -> Vec<u8> {
    let mut r = Vec::with_capacity(text.len());
    for i in 0..text.len() {
        r.push(translate_nucleotide(text[i]));
    }
    r
}

#[derive(Clone)]
pub struct Segment {
    start: usize,
    end: usize,
}
impl fmt::Debug for Segment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "S{{{}, {}}} ({})", self.start, self.end, self.end - self.start)
    }
}

fn merge_segments_with_delta(_segments: &Vec<Segment>, delta: u64) -> Vec<Segment> {
    let mut segments = _segments.clone();
    let mut r = Vec::new();
    segments.sort_by(|x, y| x.start.cmp(&y.start));

    r.push(segments[0].clone());
    for i in 1..segments.len() {
        if segments[i].start - r.last().unwrap().end <= delta as usize{
            r.last_mut().unwrap().end = segments[i].end;
        } else {
            r.push(segments[i].clone());
        }
    }
    r
}

fn merge_segments(_segments: &Vec<Segment>) -> Vec<Segment> {
    let mut segments = _segments.clone();
    let mut r = Vec::new();
    segments.sort_by(|x, y| x.start.cmp(&y.start));

    r.push(segments[0].clone());
    for i in 1..segments.len() {
        if r.last().unwrap().end < segments[i].start {
            r.push(segments[i].clone());
        } else if r.last().unwrap().end  < segments[i].end {
            r.last_mut().unwrap().end = segments[i].end;
        }
    }
    r
}

fn segments_to_segments_distance(segments: &Vec<Segment>, others: &Vec<Segment>) -> u32 {
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

pub fn search(dna: &[u8], array: &SuffixArray, pattern: &[u8], candidate_size: usize) -> Vec<Segment> {
    let mut lo = 0;
    let mut hi = dna.len()-1;
    let mut result = HashSet::new();

    while lo <= hi {
        let mid = (lo+hi)/2;
        if array[mid]+pattern.len() > dna.len() {break;}
        let substring = &dna[array[mid]..array[mid]+pattern.len()];
        let res = strcmp(&pattern, substring);

        if res > 0 {
            lo = mid + 1;
        } else if res < 0 {
            hi = mid - 1;
        } else { // TIME TO BLOB!
            let mut l = mid;
            let mut r = mid;
            let mut expand_left = true;
            let mut expand_right = true;

            while expand_left || expand_right && array[r]+pattern.len() < dna.len() {
                if expand_left {
                    if *pattern == dna[array[l]..array[l]+pattern.len()-1] {
                        result.insert(array[l]);
                        l -= 1;
                        if l == 0 { expand_left = false; }
                    } else {
                        expand_left = false;
                    }
                }
                if expand_right {
                    if *pattern == dna[array[r]..array[r]+pattern.len()] {
                        result.insert(array[r]);
                        r += 1;
                        if r == array.len() - 1 { expand_right = false; }
                    } else {
                        expand_right = false;
                    }
                }
            }
            let mut rr = Vec::new();
            for i in result {
                rr.push(Segment{start: i, end: i+candidate_size});
            }
            return rr;
        }
    }
    let mut rr = Vec::new();
    for i in result {
        rr.push(Segment{start: i, end: i+candidate_size});
    }
    return rr;
}


fn expand_nw(dna: &[u8], reverse_dna: &[u8], straight_start: usize, reverse_start: usize) -> Palindrome {
    const PRECISION: f32 = 0.8;
    const EXPANSION_STEP: usize = 1000;

    let mut offset:usize = 0;
    let mut rate = 1.0;

    let mut errors = 0;
    let mut correction_la = 0;
    let mut correction_ra = 0;

    let mut la_start;
    let mut ra_start;
    while rate >= PRECISION {
        la_start = straight_start + offset + correction_la;
        ra_start = reverse_start + offset + correction_ra;

        if (la_start+EXPANSION_STEP >= dna.len()) || (ra_start+EXPANSION_STEP >= reverse_dna.len()) { break; }

        let result = needleman_wunsch(
            &dna[la_start..la_start + EXPANSION_STEP],
            &reverse_dna[ra_start..ra_start + EXPANSION_STEP]);

        errors += result.errors;

        if result.ins_la > result.ins_ra {
            correction_ra += (result.ins_la - result.ins_ra) as usize;
        } else {
            correction_la += (result.ins_ra - result.ins_la) as usize;
        }

        rate = 1.0 - (errors as f32)/(offset as f32 + EXPANSION_STEP as f32);
        offset += EXPANSION_STEP;
    }

    Palindrome {
        left: straight_start,
        right: reverse_start,
        size: offset+EXPANSION_STEP,
        rate: rate,
    }
}
