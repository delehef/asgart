extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;
extern crate docopt;
extern crate rustc_serialize;

use std::env;
use std::io::Write;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;

use bio::io::fasta;
use bio::data_structures::suffix_array::suffix_array;

use threadpool::ThreadPool;

use docopt::Docopt;

mod utils;

const VERSION: &'static str = "
Asgart v0.1
Copyright © 2016 IRIT
";

const USAGE: &'static str = "
Usage:
    asgart <strand1-file> <strand2-file> <kmer-size> <gap-size> [-v] [-R] [-T] [-A] [--threads=<tc>] [--prefix=<prefix>]
    asgart --version
    asgart --help

Options:
    --help                  Show this message.
    --version               Show the version of Asgart.

    --threads=<tc>          Number of threads used, number of cores if 0 [default: 0].
    --prefix=<prefix>       Prefix for the result file [default: ]
    -R, --reverse           Reverse the second DNA strand.
    -T, --translate         Translate the second DNA strand.
    -A, --align             Try to perform alignment.
    -v, --verbose           Print additional informations to STDOUT.
";

#[derive(Debug, RustcDecodable)]
struct Args {
    arg_strand1_file: String,
    arg_strand2_file: String,
    arg_kmer_size: usize,
    arg_gap_size: u32,
    flag_prefix: String,
    flag_threads: usize,

    flag_reverse: bool,
    flag_translate: bool,
    flag_align: bool,

    flag_help: bool,
    flag_version: bool,
    flag_verbose: bool,
}


fn main() {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.argv(env::args()).help(true).version(Some(VERSION.to_string())).decode())
        .unwrap_or_else(|e| e.exit() )
        ;

    let out_file = format!("{}sd_{}_{}{}{}.csv",
                           args.flag_prefix,
                           args.arg_kmer_size,
                           args.arg_gap_size,
                           if args.flag_reverse {"r"} else {""},
                           if args.flag_translate {"t"} else {""},
                          );
    let threads_count: usize = if args.flag_threads > 0 { args.flag_threads } else { num_cpus::get() };

    if args.flag_verbose {
        println!("1st strand file          {}", &args.arg_strand1_file);
        println!("2nd strand file          {}", &args.arg_strand2_file);
        println!("K-mers size              {}", args.arg_kmer_size);
        println!("Max gap size             {}", args.arg_gap_size);
        println!("Output file              {}", &out_file);
        println!("Reverse 2nd strand       {}", args.flag_reverse);
        println!("Translate 2nd strand     {}", args.flag_translate);
        println!("Threads count            {}", threads_count);
        println!("");
    }

    let result = search_duplications(
        &args.arg_strand1_file, &args.arg_strand2_file,
        args.arg_kmer_size, args.arg_gap_size,
        args.flag_reverse, args.flag_translate, args.flag_align,
        threads_count
        );
    let mut out = File::create(&out_file).unwrap();
    for p in result {
        writeln!(&mut out, "{}", format!("{};{};{};{}", p.left, p.right, p.size, p.rate)).unwrap();
    }
}

fn search_duplications(
    strand1_file: &str,
    strand2_file: &str,

    kmer_size: usize,
    max_gap_size: u32,

    reverse: bool,
    translate: bool,
    align: bool,

    threads_count: usize,
    ) -> Vec<utils::SD> {

    let mut result : Vec<utils::SD> = Vec::new();

    for record1 in fasta::Reader::from_file(strand1_file).unwrap().records() {
        let strand1 = record1.unwrap().seq().to_vec();
        let shared_strand1 = Arc::new(strand1);

        for record2 in fasta::Reader::from_file(strand2_file).unwrap().records() {
            let thread_pool = ThreadPool::new(threads_count);
            let (tx, rx) = mpsc::channel();
            let strand2 = {
                let mut strand2 = record2.unwrap().seq().to_vec();
                if translate { strand2 = utils::translated(&strand2[0..strand2.len()-1].to_vec()); }
                if reverse { strand2.reverse(); }
                strand2.push(b'$');
                strand2
            };

            let shared_suffix_array = Arc::new(suffix_array(&strand2));
            let shared_strand2 = Arc::new(strand2);


            {
                const CHUNK_SIZE: usize = 1000000;
                let num_tasks = (shared_strand1.len()-kmer_size)/CHUNK_SIZE;
                let chunk_overflow = (shared_strand1.len()-kmer_size)%CHUNK_SIZE;

                let mut start = 0;
                for id in 0..num_tasks+1 // TODO Do with Vec::chunks
                {
                    let suffix_array = shared_suffix_array.clone();
                    let strand1 = shared_strand1.clone();
                    let strand2 = shared_strand2.clone();

                    let my_tx = tx.clone();

                    thread_pool.execute(move || {
                        let end = start + if id<num_tasks {CHUNK_SIZE} else {chunk_overflow};
                        my_tx.send(utils::search_duplications(
                                &strand1, &strand2, &suffix_array,
                                start, end,
                                kmer_size, max_gap_size,
                                align)).unwrap();
                    });

                    start += CHUNK_SIZE;
                }
            }

            drop(tx);
            println!("Collecting first pass...");
            let mut passes: Vec<utils::ProcessingSD> = rx.iter().fold(Vec::new(), |mut a, b| {a.append(&mut b.clone()); a});
            println!("Done.");

            if align {
                println!("Collecting second pass...");
                passes = passes.iter()
                    .map(|b| {utils::align_perfect(b.clone())}).collect();
                println!("Done.");

                println!("Collecting third pass...");
                passes = passes.iter()
                    .map(|b| {utils::align_fuzzy(&(shared_strand1.clone()), &(shared_strand2.clone()), b.clone())}).collect();
                println!("Done.");
            }

            println!("Collecting final pass...");
            result.extend(passes.iter().filter_map(|b| {
                match *b {
                    utils::ProcessingSD::Done(ref _p) => {
                        let mut p = _p.clone();
                        if reverse {
                            p.right = (*shared_strand2.clone()).len() - p.right - p.size;
                            if p.right < p.left {
                                std::mem::swap(&mut p.left, &mut p.right);
                            }
                        }
                        Some(p)
                    }
                    utils::ProcessingSD::ForFuzzy{ .. } => {println!("FOUND A FORFUZZY"); None}
                    utils::ProcessingSD::ForSW{ .. }    => {println!("FOUND A FORSW"); None},
                    utils::ProcessingSD::Empty          => {None}
                }
            }));

            // result.append(&mut all_sds);
            println!("Done for {} & {}.", kmer_size, max_gap_size);

        }
    }
    result
}
