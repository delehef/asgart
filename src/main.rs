#![allow(dead_code)]
extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;

use std::io::Write;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;

use std::str::from_utf8;

use bio::io::fasta;
use bio::data_structures::suffix_array::suffix_array;

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
    let reader = fasta::Reader::from_file("Y_soft.fasta");
    let filename = "pals.csv";
    let threads_count: usize = 3;

    println!("Threads count            {}", threads_count);
    println!("Primer size              {}", utils::CANDIDATE_SIZE);
    println!("Filename                 {}", filename);
    println!("swap_dna                 {}", cfg!(feature="swap_dna"));
    println!("");

    for record in reader.unwrap().records() {
        let threads_pool = ThreadPool::new(threads_count);
        let (tx, rx) = mpsc::channel();
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
            for id in 0..num_tasks+1 // TODO Do with Vec::chunks
            {
                println!("Launching chunk #{}", id);
                let suffix_array = shared_suffix_array.clone();
                let dna = shared_dna.clone();
                let reverse_translate_dna = shared_reverse_translate_dna.clone();

                let my_tx = tx.clone();

                threads_pool.execute(move || {
                    println!("Starting #{}", id);
                    let end = start + if id<num_tasks {CHUNK_SIZE} else {chunk_overflow};
                    if !cfg!(feature="swap_dna") {
                        my_tx.send(utils::make_palindromes(
                                &dna, &reverse_translate_dna, &suffix_array,
                                start, end)).unwrap();
                    } else {
                        my_tx.send(utils::make_palindromes(
                                &reverse_translate_dna, &dna, &suffix_array,
                                start, end)).unwrap();
                    }
                    println!("Sent for #{}", id);
                });

                start += CHUNK_SIZE;
            }
        }

        drop(tx);
        println!("Collecting first pass...");
        let first_pass: Vec<utils::ProcessingPalindrome> = rx.iter().fold(Vec::new(), |mut a, b| {a.append(&mut b.clone()); a});
        println!("Done.");

        println!("Collecting second pass...");
        let second_pass: Vec<utils::ProcessingPalindrome> = first_pass.iter()
            .map(|b| {utils::sw_align(b.clone())}).collect();
        println!("Done.");

        println!("Collecting third pass...");
        let third_pass: Vec<utils::ProcessingPalindrome> = second_pass.iter()
            .map(|b| {utils::fuzzy_align(&(shared_dna.clone()), &(shared_reverse_translate_dna.clone()), b.clone())}).collect();
        println!("Done.");

        println!("Collecting final pass...");
        let final_pass = third_pass.iter()
            .map(|b| {
                match *b {
                    utils::ProcessingPalindrome::Done(ref p) => Some(p),
                    utils::ProcessingPalindrome::ForFuzzy{pp:_, right:_} => {println!("FOUND A FORFUZZY"); None}
                    utils::ProcessingPalindrome::ForSW{pp: _, left_match: _, right_match: _, right_segment: _} => {println!("FOUND A FORSW"); None},
                    utils::ProcessingPalindrome::TooSmall => {None}
                    utils::ProcessingPalindrome::Empty => {None}
                }
            });
        println!("Done.");

        let mut out = File::create(filename).unwrap();
        for palindrome in final_pass {
            let to_write = palindrome.map_or((), |p| { writeln!(&mut out, "{}", format!("{};{};{};{}", p.left, p.right, p.size, p.rate)).ok().unwrap()});
        }
    }
}
