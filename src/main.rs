#![feature(iter_cmp)]
extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;

use std::io::prelude::*;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;

use std::collections::HashSet;
use std::str::from_utf8;
use std::cmp;

use bio::io::fasta;
use bio::data_structures::suffix_array::SuffixArray;
use bio::data_structures::suffix_array::suffix_array;

use uuid::Uuid;

use threadpool::ThreadPool;

mod utils;


const CANDIDATE_SIZE: usize = 20;
const PALINDROME_THRESHOLD_SIZE: usize = 1000;
const PRIMER_SHIFT: usize = 10;

#[derive(Debug)]
struct Palindrome {
    left: usize,
    right: usize,
    size: usize,
    rate: f32,
}



fn expand_nw(dna: &[u8], reverse_dna: &[u8], straight_start: usize, reverse_start: usize) -> Palindrome {
    const PRECISION: f32 = 0.9;
    const EXPANSION_STEP: usize = 1000;

    let mut offset:usize = 0;
    let mut rate = 1.0;

    let mut errors = 0;
    let mut correction_la = 0;
    let mut correction_ra = 0;

    let mut la_start;
    let mut ra_start;
    while rate >= PRECISION {
        offset += EXPANSION_STEP;
        la_start = straight_start + offset + correction_la;
        ra_start = reverse_start + offset + correction_ra;

        if (la_start+EXPANSION_STEP >= dna.len()) || (ra_start+EXPANSION_STEP >= reverse_dna.len()) { break; }

        let result = utils::needleman_wunsch(
            &dna[la_start..la_start + EXPANSION_STEP],
            &reverse_dna[ra_start..ra_start + EXPANSION_STEP]);

        errors += result.errors;

        if result.ins_la > result.ins_ra {
            correction_ra += (result.ins_la - result.ins_ra) as usize;
        } else {
            correction_la += (result.ins_ra - result.ins_la) as usize;
        }

        rate = 1.0 - (errors as f32)/(offset as f32 + EXPANSION_STEP as f32);
    }

    println!("{};{};{};{}",
             straight_start,
             reverse_start,
             offset+EXPANSION_STEP,
             rate
            );

    Palindrome {
        left: straight_start,
        right: reverse_start,
        size: offset+EXPANSION_STEP,
        rate: rate,
    }
}

fn expand_palindrome(dna: &[u8], reverse_dna: &[u8], straight_start: usize, reverse_start: usize) -> Palindrome {
    const PRECISION: f32 = 0.9;
    const EXPANSION_STEP: usize = 100;

    let mut expansion = 1;
    let mut current_rate = 1.0;

    while current_rate > PRECISION && straight_start+expansion < dna.len() - EXPANSION_STEP {
        expansion += EXPANSION_STEP;
        let mut mutations = 0;
        for i in 0..expansion {
            let na = dna[straight_start + i];
            let nb = reverse_dna[reverse_start + i];

            if na != nb { mutations += 1; }
        }
        current_rate = 1.0 - (mutations as f32/expansion as f32);
    }

    let straight_candidate: usize = dna.len() - (reverse_start-CANDIDATE_SIZE);

    Palindrome {
        left: cmp::min(straight_start, straight_candidate),
        right: cmp::max(straight_start, straight_candidate),
        size: expansion,
        rate: current_rate,
    }
}

fn search(dna: &[u8], array: &SuffixArray, pattern: &[u8]) -> HashSet<usize> {
    let mut lo = 0;
    let mut hi = dna.len()-1;
    let mut result = HashSet::new();

    while lo <= hi {
        let mid = (lo+hi)/2;
        if array[mid]+pattern.len() > dna.len() {break;}
        let substring = &dna[array[mid]..array[mid]+pattern.len()];
        let res = utils::strcmp(&pattern, substring);

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

fn window(text: &[u8], begin: usize, extension: usize) -> &str {
    if begin+extension <= text.len() {
        from_utf8(&text[begin..begin+extension]).unwrap()
    } else {
        "N/A"
    }
}

fn look_for_palindromes(dna: &[u8], reverse_translate_dna: &[u8], sa: &SuffixArray, start: usize, end: usize) -> Vec<Palindrome> {
    let mut palindromes = Vec::new();

    let mut i = start;
    while i < end {
        let bottom = i;
        let top = bottom + CANDIDATE_SIZE;
        let candidate = &dna[bottom..top];
        if candidate[0] == 'N' as u8 { continue;  }
        let results = search(&reverse_translate_dna, &sa, candidate);
        if results.len() > 0 {
            //let mut min_size = dna.len();
            for result in results {
                let mut palindrome = expand_palindrome(&dna, &reverse_translate_dna, bottom, result);
                if palindrome.size >= PALINDROME_THRESHOLD_SIZE {
                    palindrome = expand_nw(&dna, &reverse_translate_dna, bottom, result);
                    //if palindrome.size < min_size { min_size = palindrome.size; }
                    palindromes.push(palindrome);
                }
            }
            //if min_size == dna.len() { min_size = 0 };
            //i += min_size;
        }
        i += PRIMER_SHIFT;
    }

    palindromes
}

fn main () {
    let reader = fasta::Reader::from_file("Y.fasta");
    let filename = format!("pals-{}.csv", Uuid::new_v4().to_string());
    let threads_count: usize = 3; //num_cpus::get();
    let threads_pool = ThreadPool::new(threads_count);
    println!("Threads count            {}", threads_count);
    println!("Primer size              {}", CANDIDATE_SIZE);
    println!("NW extensions threshold  {}", PALINDROME_THRESHOLD_SIZE);
    println!("Primer window base shift {}", PRIMER_SHIFT);
    println!("Filename                 {}", filename);
    println!("swap_dna                 {}", cfg!(feature="swap_dna"));
    println!("");

    let mut out = File::create(filename).unwrap();
    let (tx, rx) = mpsc::channel();
    for record in reader.unwrap().records() {
        let seqs = record.unwrap();
        let dna = seqs.seq().to_vec();
        let shared_dna = Arc::new(dna);
        let mut reverse_translate_dna = utils::translate(&shared_dna[0..shared_dna.len()-1].to_vec());
        reverse_translate_dna.reverse();
        reverse_translate_dna.push('$' as u8);

        let shared_sa = Arc::new(suffix_array(&reverse_translate_dna));
        let shared_reverse_translate_dna = Arc::new(reverse_translate_dna);


        {
            const CHUNK_SIZE: usize = 200000;
            let num_tasks = (shared_dna.len()-CANDIDATE_SIZE)/CHUNK_SIZE;
            let chunk_overflow = (shared_dna.len()-CANDIDATE_SIZE)%CHUNK_SIZE;

            let mut start = 19000000;
            for id in 0..num_tasks+1
            {
                let local_sa = shared_sa.clone();
                let local_dna = shared_dna.clone();
                let local_reverse_translate_dna = shared_reverse_translate_dna.clone();

                let my_tx = tx.clone();

                threads_pool.execute(move || {
                    let end = start + if id<num_tasks {CHUNK_SIZE} else {chunk_overflow};
                    // println!("Going from {} to {}", start, end);
                    // my_tx.send(look_for_palindromes(
                    //         &local_dna, &local_reverse_translate_dna, &local_sa,
                    //         start, end)).unwrap();
                    if !cfg!(feature="swap_dna") {
                        my_tx.send(make_palindromes(
                                &local_dna, &local_reverse_translate_dna, &local_sa,
                                start, end)).unwrap();
                    } else {
                        my_tx.send(make_palindromes(
                                &local_reverse_translate_dna, &local_dna, &local_sa,
                                start, end)).unwrap();
                    }
                });

                start += CHUNK_SIZE;
            }
        }
    }

    drop(tx);
    for palindromes_list in rx.iter() {
        for palindrome in palindromes_list {
            writeln!(&mut out,
                     "{};{};{};{}",
                     palindrome.left,
                     palindrome.right,
                     palindrome.size,
                     palindrome.rate
                    ).unwrap();
        }
    }
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

const M: u32 = 3*CANDIDATE_SIZE as u32;
const MAX_HOLE_SIZE: u32 = 1000;
struct ProtoPalindrome {
    bottom: usize,
    top: usize,
    matches: Vec<HashSet<usize>>,
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


    let (right_begin, _) = *right_arms.iter().max_by(|&&(x, y)| y-x).unwrap();

    Palindrome {
        left: pp.bottom,
        right: *right_begin,
        size: pp.top - pp.bottom,
        rate: 0.0
    }
}

enum SearchState {
    START,
    Grow,
    SparseGrow,
    PROTO,
}

fn make_palindromes(dna: &[u8], rt_dna: &[u8], sa: &SuffixArray, start: usize, end: usize) -> Vec<Palindrome> {
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
                    // println!("Starting @{}", i);
                    current_sets.push(new_set);
                    state = SearchState::Grow;
                } else {
                    // println!("Nothing @{}", i);
                    i += 1;
                    state = SearchState::START;
                }

            },
            SearchState::Grow => {
                // println!("Growing @{}", i);
                i += CANDIDATE_SIZE/5;
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
                // println!("Sparse @{}", i);
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

