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




fn window(text: &[u8], begin: usize, extension: usize) -> &str {
    if begin+extension <= text.len() {
        from_utf8(&text[begin..begin+extension]).unwrap()
    } else {
        "N/A"
    }
}



fn main () {
    let reader = fasta::Reader::from_file("Y.fasta");
    let filename = "pals.csv";
    let threads_count: usize = 3;
    let threads_pool = ThreadPool::new(threads_count);

    println!("Threads count            {}", threads_count);
    println!("Primer size              {}", utils::CANDIDATE_SIZE);
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
