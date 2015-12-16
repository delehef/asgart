use std::collections::HashSet;
use std::cmp;
use bio::data_structures::suffix_array::SuffixArray;
use std::io::Write;

pub const M: u32 = 3*CANDIDATE_SIZE as u32;
pub const MAX_HOLE_SIZE: u32 = 1000;
pub const CANDIDATE_SIZE: usize = 20;


macro_rules! log(
    ($($arg:tt)*) => (
        match writeln!(&mut ::std::io::stderr(), $($arg)* ) {
            Ok(_) => {},
            Err(x) => panic!("Unable to write to stderr: {}", x),
        }
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
    START,
    Grow,
    SparseGrow,
    PROTO,
}

struct ProtoPalindrome {
    bottom: usize,
    top: usize,
    matches: Vec<HashSet<usize>>,
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

fn make_palindrome(pp: &ProtoPalindrome) -> Palindrome {
    let mut matches = Vec::new();
    for set in pp.matches.iter() {
        for k in set {
            matches.push(k);
        }
    }
    matches.sort();
    matches.dedup();

    let mut right_arms = Vec::new();
    let mut current_start = matches[0];
    for i in 1..matches.len() {
        // println!("{}", matches[i]);
        if matches[i] - matches[i-1] > 20000 as usize {
            right_arms.push((current_start, matches[i-1]+CANDIDATE_SIZE));

            current_start = matches[i];
        }
    }
    // print!("Possibilities: {:?}", right_arms);


    let (right_begin, _) = *right_arms.iter().max_by_key(|&&(x, y)| y-x).unwrap();

    Palindrome {
        left: pp.bottom,
        right: *right_begin,
        size: pp.top - pp.bottom,
        rate: 0.0
    }
}

pub fn make_palindromes(dna: &[u8], rt_dna: &[u8], sa: &SuffixArray, start: usize, end: usize) -> Vec<Palindrome> {
    let mut r = Vec::new();

    let mut i = start;
    let mut hole = 0;
    let mut current_start = 0;
    let mut current_sets = Vec::new();
    let mut state = SearchState::START;

    loop {
        match state {
            SearchState::START => {
                if i+1 >= end {
                    break;
                }
                hole = 0;
                current_start = i;
                current_sets = Vec::new();
                let new_set = search(&rt_dna, &sa, &dna[i..i+CANDIDATE_SIZE]);
                if new_set.len() > 0 {
                    log!("Starting @{}", i);
                    current_sets.push(new_set);
                    state = SearchState::Grow;
                } else {
                    log!("Nothing @{}", i);
                    i += 1;
                    state = SearchState::START;
                }

            },
            SearchState::Grow => {
                log!("Growing @{}", i);
                i += CANDIDATE_SIZE/4;
                let set = search(&rt_dna, &sa, &dna[i..i+CANDIDATE_SIZE]);

                if i >= dna.len() - CANDIDATE_SIZE {
                    state = SearchState::PROTO;
                } else if set_to_sets_distance(&set, &current_sets) <= 20000/*M*hole*/ {
                    current_sets.push(set);
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::SparseGrow => {
                log!("Sparse @{}", i);
                i += 1;
                hole += 1;
                let set = search(&rt_dna, &sa, &dna[i..i+CANDIDATE_SIZE]);

                if hole > MAX_HOLE_SIZE {
                    state = SearchState::PROTO;
                } else if i >= dna.len() - CANDIDATE_SIZE {
                    state = SearchState::PROTO;
                } else if set_to_sets_distance(&set, &current_sets) <= 20000/*M*hole*/ {
                    current_sets.push(set);
                    state = SearchState::Grow;
                } else {
                    state = SearchState::SparseGrow;
                }
            },
            SearchState::PROTO => {
                if i - current_start > 10000 {
                    let pp = ProtoPalindrome {
                        bottom: current_start,
                        top: i,
                        matches: current_sets.clone()
                    };
                    let p = make_palindrome(&pp);
                    println!("{};{};{}", p.left, p.right, p.size);
                    r.push(p);
                }
                state = SearchState::START;
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


pub fn search(dna: &[u8], array: &SuffixArray, pattern: &[u8]) -> HashSet<usize> {
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
            return result;
        }
    }
    result
}
