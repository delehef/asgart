#![feature(iter_cmp)]
#![allow(dead_code)]
extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;

use std::io::prelude::*;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;

use std::str::from_utf8;
use std::cmp;

use bio::io::fasta;
use bio::data_structures::suffix_array::suffix_array;

use uuid::Uuid;

use threadpool::ThreadPool;

mod utils;


const PALINDROME_THRESHOLD_SIZE: usize = 1000;
const PRIMER_SHIFT: usize = 10;




fn expand_nw(dna: &[u8], reverse_dna: &[u8], straight_start: usize, reverse_start: usize) -> utils::Palindrome {
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

    utils::Palindrome {
        left: straight_start,
        right: reverse_start,
        size: offset+EXPANSION_STEP,
        rate: rate,
    }
}

fn expand_palindrome(dna: &[u8], reverse_dna: &[u8], straight_start: usize, reverse_start: usize) -> utils::Palindrome {
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

    let straight_candidate: usize = dna.len() - (reverse_start-utils::CANDIDATE_SIZE);

    utils::Palindrome {
        left: cmp::min(straight_start, straight_candidate),
        right: cmp::max(straight_start, straight_candidate),
        size: expansion,
        rate: current_rate,
    }
}

fn window(text: &[u8], begin: usize, extension: usize) -> &str {
    if begin+extension <= text.len() {
        from_utf8(&text[begin..begin+extension]).unwrap()
    } else {
        "N/A"
    }
}

//fn look_for_palindromes(dna: &[u8], reverse_translate_dna: &[u8], sa: &SuffixArray, start: usize, end: usize) -> Vec<utils::Palindrome> {
//    let mut palindromes = Vec::new();

//    let mut i = start;
//    while i < end {
//        let bottom = i;
//        let top = bottom + utils::CANDIDATE_SIZE;
//        let candidate = &dna[bottom..top];
//        if candidate[0] == 'N' as u8 { continue;  }
//        let results = utils::search(&reverse_translate_dna, &sa, candidate);
//        if results.len() > 0 {
//            //let mut min_size = dna.len();
//            for result in results {
//                let mut palindrome = expand_palindrome(&dna, &reverse_translate_dna, bottom, result);
//                if palindrome.size >= PALINDROME_THRESHOLD_SIZE {
//                    palindrome = expand_nw(&dna, &reverse_translate_dna, bottom, result);
//                    //if palindrome.size < min_size { min_size = palindrome.size; }
//                    palindromes.push(palindrome);
//                }
//            }
//            //if min_size == dna.len() { min_size = 0 };
//            //i += min_size;
//        }
//        i += PRIMER_SHIFT;
//    }

//    palindromes
//}


fn main () {
    let reader = fasta::Reader::from_file("Y_soft.fasta");
    let filename = format!("pals-{}.csv", Uuid::new_v4().to_string());
    let threads_count: usize = num_cpus::get()/2;
    let threads_pool = ThreadPool::new(threads_count);

    println!("Threads count            {}", threads_count);
    println!("Primer size              {}", utils::CANDIDATE_SIZE);
    println!("NW extensions threshold  {}", PALINDROME_THRESHOLD_SIZE);
    println!("Primer window base shift {}", PRIMER_SHIFT);
    println!("Filename                 {}", filename);
    println!("swap_dna                 {}", cfg!(feature="swap_dna"));
    println!("");

    let (tx, rx) = mpsc::channel();
    for record in reader.unwrap().records() {
        let dna = record.unwrap().seq().to_vec();
        // Build the RT DNA sequence
        let mut reverse_translate_dna = utils::translated(&dna[0..dna.len()-1].to_vec());
        reverse_translate_dna.reverse();
        reverse_translate_dna.push('$' as u8);

        // Shared RO suffix array
        let shared_suffix_array = Arc::new(suffix_array(&reverse_translate_dna));
        // Sharerd RO DNA sequence
        let shared_dna = Arc::new(dna);
        // Final RO reverse-translated DNA sequence
        let shared_reverse_translate_dna = Arc::new(reverse_translate_dna);


        {
            const CHUNK_SIZE: usize = 200000;
            let num_tasks = (shared_dna.len()-utils::CANDIDATE_SIZE)/CHUNK_SIZE;
            let chunk_overflow = (shared_dna.len()-utils::CANDIDATE_SIZE)%CHUNK_SIZE;

            let mut start = 0;
            for id in 0..num_tasks+1
            {
                let suffix_array = shared_suffix_array.clone();
                let dna = shared_dna.clone();
                let reverse_translate_dna = shared_reverse_translate_dna.clone();

                let my_tx = tx.clone();

                threads_pool.execute(move || {
                    let end = start + if id<num_tasks {CHUNK_SIZE} else {chunk_overflow};
                    // my_tx.send(look_for_palindromes(
                    //         &local_dna, &local_reverse_translate_dna, &local_sa,
                    //         start, end)).unwrap();
                    if !cfg!(feature="swap_dna") {
                        my_tx.send(utils::make_palindromes(
                                &dna, &reverse_translate_dna, &suffix_array,
                                start, end)).unwrap();
                    } else {
                        my_tx.send(utils::make_palindromes(
                                &reverse_translate_dna, &dna, &suffix_array,
                                start, end)).unwrap();
                    }
                });

                start += CHUNK_SIZE;
            }
        }
    }

    drop(tx);
    let mut out = File::create(filename).unwrap();
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
