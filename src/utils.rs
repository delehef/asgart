use std::collections::HashSet;
use std::cmp;
use bio::data_structures::suffix_array::SuffixArray;
use std::io::Write;
use std::fmt;
use bio::alignment::pairwise::*;

pub const M: u32 = 3*CANDIDATE_SIZE as u32;
pub const MAX_HOLE_SIZE: u32 = 2000;
pub const CANDIDATE_SIZE: usize = 20;
const MIN_PALINDROME_SIZE: usize = 1000;


macro_rules! log(
    ($($arg:tt)*) => (
        // match writeln!(&mut ::std::io::stderr(), $($arg)* ) {
        //     Ok(_) => {},
        //     Err(x) => {} //panic!("Unable to write to stderr: {}", x),
        // }
    )
);


#[derive(Debug)]
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

struct ProtoPalindrome {
    bottom: usize,
    top: usize,
    matches: Vec<Segment>,
}

enum ProcessingPalindrome {
    Done(Palindrome),
    TooLong{left: usize, right: usize, size: usize},
    TooSmall,
    Empty,
}

fn set_to_sets_distance(s: &HashSet<usize>, t: &Vec<HashSet<usize>>) -> u32 {
    let mut min = 100000000;
    for i in s {
        for ss in t {
            for j in ss {
                let diff = (*j as i32 - *i as i32).abs() as u32;
                if diff < min { min = diff }
            }
        }
    }
    min
}

fn make_palindrome(pp: &ProtoPalindrome, dna: &[u8], rt_dna: &[u8]) -> ProcessingPalindrome {
    let mut matches = pp.matches.clone();
    matches = merge_segments_with_delta(&matches, MAX_HOLE_SIZE as u64);

    // Sort by size, descending
    matches.sort_by(|a, b| (b.end - b.start).cmp(&(a.end - a.start)));
    if matches[0].end - matches[0].start < MIN_PALINDROME_SIZE { return ProcessingPalindrome::TooSmall; }

    // Fetch left and right candidate areas
    let left_match = &dna[pp.bottom..pp.top];
    let right_match = &rt_dna[matches[0].start..matches[0].end];

    if right_match.len() > 75000 || left_match.len() > 75000 {
        return ProcessingPalindrome::TooLong {
            left: pp.bottom,
            right: matches[0].start,
            size: cmp::max(right_match.len(), left_match.len())
        };
    }

    // Align them if not too big
    let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
    let mut aligner = Aligner::with_capacity(left_match.len(), right_match.len(), -5, -1, &score);
    let alignment = aligner.local(&left_match, &right_match);
    if alignment.operations.len() < MIN_PALINDROME_SIZE {return ProcessingPalindrome::Empty;};

    let result = Palindrome {
        left: pp.bottom + alignment.xstart,
        right: matches[0].start + alignment.ystart,
        size: alignment.operations.len(),
        rate: 0.0
    };
    println!("Left position:  {}", result.left);
    println!("Right position: {}", result.right);
    println!("Size:           {}", result.size);
    println!("Alignment: \n{}\n", alignment.pretty(left_match, right_match));

    return ProcessingPalindrome::Done(result);
}

pub fn make_palindromes(dna: &[u8], rt_dna: &[u8], sa: &SuffixArray, start: usize, end: usize) -> Vec<Palindrome> {
    let mut r = Vec::new();

    let mut i = start;
    let mut hole = 0;
    let mut current_start = 0;
    let mut current_segments = Vec::new();
    let mut state = SearchState::Start;

    loop {
        match state {
            SearchState::Start => {
                if i+1 >= end {
                    break;
                }
                hole = 0;
                current_start = i;
                current_segments = Vec::new();
                let mut new_segments = search(&rt_dna, &sa, &dna[i..i+CANDIDATE_SIZE]);
                if new_segments.len() > 0 {
                    log!("Starting @{}", i);
                    current_segments = merge_segments(&mut new_segments);
                    state = SearchState::Grow;
                } else {
                    i += 1;
                    state = SearchState::Start;
                }

            },
            SearchState::Grow => {
                log!("Growing @{}", i);
                i += CANDIDATE_SIZE/2;
                let mut new_matches = search(&rt_dna, &sa, &dna[i..i+CANDIDATE_SIZE]);

                if i >= dna.len() - CANDIDATE_SIZE {
                    state = SearchState::Proto;
                } else if segments_to_segments_distance(&new_matches, &current_segments) <= MAX_HOLE_SIZE {
                    current_segments.append(&mut new_matches);
                    current_segments = merge_segments(&mut current_segments);
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::SparseGrow => {
                log!("Sparse @{}", i);
                i += 1;
                hole += 1;
                let mut new_matches = search(&rt_dna, &sa, &dna[i..i+CANDIDATE_SIZE]);

                if hole > MAX_HOLE_SIZE {
                    state = SearchState::Proto;
                } else if i >= dna.len() - CANDIDATE_SIZE {
                    state = SearchState::Proto;
                } else if segments_to_segments_distance(&new_matches, &current_segments) <= MAX_HOLE_SIZE {
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
                    match make_palindrome(&pp, dna, rt_dna) {
                        ProcessingPalindrome::Done(p) => {
                            println!("{};{};{}", p.left, p.right, p.size);
                            r.push(p);
                        },
                        ProcessingPalindrome::TooLong{left, right, size} => {
                            println!("Too long @{}-{}/{}", left, right, size);
                        }
                        _ => {},
                    }
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

pub fn search(dna: &[u8], array: &SuffixArray, pattern: &[u8]) -> Vec<Segment> {
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
                rr.push(Segment{start: i, end: i+CANDIDATE_SIZE});
            }
            return rr;
        }
    }
    let mut rr = Vec::new();
    for i in result {
        rr.push(Segment{start: i, end: i+CANDIDATE_SIZE});
    }
    return rr;
}
