#![feature(unicode)]
extern crate bio;
use std::io::prelude::*;
use std::fs::File;

use std::collections::HashSet;
use std::str::from_utf8;
use std::cmp;

use std::cmp::Ordering;
use bio::io::fasta;
use bio::data_structures::suffix_array::SuffixArray;
use bio::data_structures::suffix_array::suffix_array;


const CANDIDATE_SIZE: usize = 20;

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

fn needleman_wunsch(q1: &[u8], q2: &[u8]) -> i32 {
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

    let mut alignment_a = String::with_capacity(q1.len());
    let mut alignment_b = String::with_capacity(q2.len());
    let mut alignment_f = String::with_capacity(q2.len());
    let mut error = 0;

    let mut result = AlignmentResult {errors: 0, ins_la: 0, ins_ra: 0};
    while i>1 && j>1 {
        if i>0 && j>0
            && (F[i+j*size] == F[(i-1)+(j-1)*size] + if q1[i] == q2[j] {MATCH_SCORE} else {MISMATCH_SCORE})
        {
            // alignment_a.push(q1[i] as char);
            // alignment_b.push(q2[j] as char);
            if q1[i] != q2[j] {
                result.errors += 1;
                // alignment_f.push('#');
            } else {
                // alignment_f.push(' ');
            }
            i -= 1;
            j -= 1;
        } else if i>0 && F[i+j*size] == F[(i-1)+j*size] + DELETE_SCORE {
            // alignment_a.push(q1[i] as char);
            // alignment_b.push('-');
            // alignment_f.push('-');
            result.ins_la += 1;
            i -= 1;
        } else if j>0 && F[i+j*size] == F[i+(j-1)*size] + DELETE_SCORE {
            // alignment_a.push('-');
            // alignment_b.push(q2[j] as char);
            // alignment_f.push('-');
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
    const precision: f32 = 0.9;
    const expansion_step: usize = 1000;
    const window_size: usize = 3000;

    let mut expansion:usize = 0;
    let mut rate = 1.0;

    while rate >= precision {
        expansion += expansion_step;
        let result = needleman_wunsch(
            &dna[straight_start..straight_start+expansion],
            &reverse_dna[reverse_start..reverse_start+expansion]);
        rate = 1.0 - (errors as f32)/(expansion as f32);
    }

    Palindrome {
        left: straight_start,
        right: reverse_start,
        size: expansion,
        rate: rate,
    }
}

fn expand_palindrome(dna: &[u8], reverse_dna: &[u8], straight_start: usize, reverse_start: usize) -> Palindrome {
    const precision: f32 = 0.9;
    const expansion_step: usize = 100;

    let mut expansion = 1;
    let mut current_rate = 1.0;

    while current_rate > precision && straight_start+expansion < dna.len() - expansion_step {
        expansion += expansion_step;
        let mut sa = Vec::new();
        let mut sb = Vec::new();

        let mut mutations = 0;
        for i in 0..expansion {
            let na = dna[straight_start + i];
            let nb = reverse_dna[reverse_start + i];
            sa.push(na);
            sb.push(nb);

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
        //println!("\n\n{} {} {}", lo, mid, hi);
        if array[mid]+pattern.len() > dna.len() {break;}
        let substring = &dna[array[mid]..array[mid]+pattern.len()];
        let res = strcmp(&pattern, substring);
        //println!("res = {}", res);

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

fn main ()
{
    // let mut f1 = File::open("/home/franklin/bioinfo/needleman_wunsch/test1_la.fasta").unwrap();
    // let mut f2 = File::open("/home/franklin/bioinfo/needleman_wunsch/test1_ra.fasta").unwrap();
    // let mut la = String::new();
    // let mut ra = String::new();
    // f1.read_to_string(&mut la);
    // f2.read_to_string(&mut ra);
    // println!("{}", needleman_wunsch(la.as_bytes(), ra.as_bytes()));

    let reader = fasta::Reader::from_file("/home/franklin/Y.fasta");

    for record in reader.unwrap().records() {
        let seqs = record.unwrap();
        let dna = seqs.seq();
        let mut reverse_translate_dna = translate(&dna[0..dna.len()-1].to_vec());
        reverse_translate_dna.reverse();
        reverse_translate_dna.push('$' as u8);

        let sa = suffix_array(&reverse_translate_dna);

        let mut i = 0;
        while i < dna.len()-CANDIDATE_SIZE {
            let bottom = i;
            let top = bottom + CANDIDATE_SIZE;
            let candidate = &dna[bottom..top];
            if candidate == "CATTGTAGTTAATGCCCTGA".as_bytes() { println!("Looking for P6"); }
            let results = search(&reverse_translate_dna, &sa, candidate);
            if results.len() > 0 {
                // println!("{} -> {}", bottom, top);
                // println!("{:?}", results);
                let mut min_size = 0;
                for result in results {
                    let mut palindrome = expand_palindrome(dna, &reverse_translate_dna, bottom, result);
                    if palindrome.size > 1000 {
                        palindrome = expand_nw(dna, &reverse_translate_dna, bottom, result);
                        println!("{};{};{};{}", palindrome.left, palindrome.right, palindrome.size, palindrome.rate);

                        if min_size == 0 { min_size = palindrome.size }
                        else if palindrome.size < min_size { min_size = palindrome.size; }
                    }
                }
                i += min_size;
            }
            i += 1;
        }
    }
}
