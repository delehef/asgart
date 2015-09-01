extern crate uuid;
extern crate bio;
extern crate num_cpus;

use std::io::prelude::*;
use std::fs::File;
use std::thread;
use std::sync::mpsc;
use std::sync::Arc;

use std::collections::HashSet;
use std::str::from_utf8;
use std::cmp;

use bio::io::fasta;
use bio::data_structures::suffix_array::SuffixArray;
use bio::data_structures::suffix_array::suffix_array;

use uuid::Uuid;


const CANDIDATE_SIZE: usize = 20;
const PALINDROME_THRESHOLD_SIZE: usize = 1000;
const PRIMER_SHIFT: usize = 1;

#[derive(Debug)]
struct Palindrome {
    left: usize,
    right: usize,
    size: usize,
    rate: f32,
}

struct AlignmentResult {
    errors: u32,
    ins_la: u32,
    ins_ra: u32,
}

fn needleman_wunsch(q1: &[u8], q2: &[u8]) -> AlignmentResult {
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

fn strcmp(s1: &[u8], s2: &[u8]) -> i32 {
    assert!(s1.len() == s2.len());
    let n = s2.len();

    for i in 0..n {
        if s1[i] != s2[i] {
            return s1[i] as i32 - s2[i] as i32;
        }
    }
    0
}

fn translate_nucleotide(n: u8) -> u8 {
    if n == 'A' as u8 { 'T' as u8}
    else if n == 'T' as u8 { 'A' as u8 }
    else if n == 'G' as u8 { 'C' as u8 }
    else if n == 'C' as u8 { 'G' as u8 }
    else { panic!("Not DNA: {}", n as u8); }
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

fn translate(text: &[u8]) -> Vec<u8> {
    let mut r = Vec::with_capacity(text.len());
    for i in 0..text.len() {
        r.push(translate_nucleotide(text[i]));
    }
    r
}

fn search(dna: &[u8], array: &SuffixArray, pattern: &[u8]) -> HashSet<usize> {
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

fn main ()
{
    let reader = fasta::Reader::from_file("Y.fasta");
    let filename = format!("pals-{}.csv", Uuid::new_v4().to_string());
    let threads_count: usize = num_cpus::get();
    println!("Threads count            {}", threads_count);
    println!("Primer size              {}", CANDIDATE_SIZE);
    println!("NW extensions threshold  {}", PALINDROME_THRESHOLD_SIZE);
    println!("Primer window base shift {}", PRIMER_SHIFT);
    println!("Filename                 {}\n", filename);

    let mut out = File::create(filename).unwrap();
    let (tx, rx) = mpsc::channel();
    for record in reader.unwrap().records() {
        let seqs = record.unwrap();
        let dna = seqs.seq().to_vec();
        let shared_dna = Arc::new(dna);
        let mut reverse_translate_dna = translate(&shared_dna[0..shared_dna.len()-1].to_vec());
        reverse_translate_dna.reverse();
        reverse_translate_dna.push('$' as u8);

        let shared_sa = Arc::new(suffix_array(&reverse_translate_dna));
        let shared_reverse_translate_dna = Arc::new(reverse_translate_dna);


        {
            let num_tasks_per_thread = (shared_dna.len()-CANDIDATE_SIZE) / threads_count;
            let num_tougher_threads = (shared_dna.len()-CANDIDATE_SIZE) % threads_count;
            let mut start = 0;
            for id in 0..threads_count {
                let local_sa = shared_sa.clone();
                let local_dna = shared_dna.clone();
                let local_reverse_translate_dna = shared_reverse_translate_dna.clone();

                let chunksize =
                    if id < num_tougher_threads {
                        num_tasks_per_thread + 1
                    } else {
                        num_tasks_per_thread
                    };
                let my_tx = tx.clone();
                thread::spawn(move || {
                    let end = start+chunksize;
                    println!("Going from {} to {}", start, end);
                    my_tx.send(look_for_palindromes(
                            &local_dna, &local_reverse_translate_dna, &local_sa,
                            start, end)).unwrap();
                });
                start += chunksize;
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
