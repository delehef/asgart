extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;
extern crate docopt;
extern crate rustc_serialize;

use std::cmp;
use std::env;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;
use std::ascii::AsciiExt;

use threadpool::ThreadPool;

use docopt::Docopt;

use divsufsort64::idx;

mod utils;
mod divsufsort64;

const VERSION: &'static str = "
Asgart v0.1
Copyright © 2016 IRIT
";

const USAGE: &'static str = "
Usage:
    asgart <strand1-file> <strand2-file> <kmer-size> <gap-size> [-v] [-R] [-T] [-A] [-i] [--threads=<tc>] [--prefix=<prefix>]
    asgart --version
    asgart --help

Options:
    --help                  Show this message.
    --version               Show the version of Asgart.

    --threads=<tc>          Number of threads used, number of cores if 0 [default: 0].
    --prefix=<prefix>       Prefix for the result file [default: ]
    -R, --reverse           Reverse the second DNA strand.
    -T, --translate         Translate the second DNA strand.
    -i, --interlaced        Look for interlaced duplications (may drastically impair performances!)
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
    flag_interlaced: bool,

    flag_help: bool,
    flag_version: bool,
    flag_verbose: bool,
}


fn read_fasta(filename: &str) -> Result<Vec<u8>, io::Error> {
    let file = try!(File::open(filename));
    let file = io::BufReader::new(file);


    let mut r = file.lines()
       .filter_map(|result| result.ok())
       .filter(|line| !line.starts_with('>'))
       .fold(Vec::new(), |mut r, line| {r.extend(line.trim().as_bytes().iter().cloned()); r});

    r = r.to_ascii_uppercase();

    // r.retain(|c| *c != b'n' && *c != b'N');

    Ok(r)
}

pub fn r_divsufsort(dna: &[u8]) -> Vec<idx> {
    let mut sa = Vec::with_capacity(dna.len());
    sa.resize(dna.len(), 0);
    unsafe {
        divsufsort64::divsufsort64(dna.as_ptr(), sa.as_mut_ptr(), dna.len() as i64);
    }
    sa
}


#[derive(RustcEncodable)]
struct Strand {
    length: usize,
    reversed: bool,
    translated: bool,
}

#[derive(RustcEncodable)]
struct RunResult {
    strand1: Strand,
    strand2: Strand,

    SDs: Vec<utils::SD>,
}

fn main() {
    let args: Args = Docopt::new(USAGE)
        .and_then(|d| d.argv(env::args()).help(true).version(Some(VERSION.to_string())).decode())
        .unwrap_or_else(|e| e.exit() )
        ;

    let out_file = format!("{}{}_vs_{}_{}_{}{}{}.json",
                           args.flag_prefix,
                           &args.arg_strand1_file,
                           &args.arg_strand2_file,
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
        println!("Interlaced SD            {}", args.flag_interlaced);
        println!("Threads count            {}", threads_count);
        println!("libdivsufsort            v{:?}", unsafe{divsufsort64::divsufsort64_version()});
        println!("");
    }

    let result = search_duplications(
        &args.arg_strand1_file, &args.arg_strand2_file,
        args.arg_kmer_size, args.arg_gap_size + args.arg_kmer_size as u32,
        args.flag_reverse, args.flag_translate, args.flag_align, args.flag_interlaced,
        threads_count
        );
    let mut out = File::create(&out_file).expect(&format!("Unable to create `{}` for output", 
                                                          &out_file));
    writeln!(&mut out, "{}", rustc_serialize::json::encode(&result).unwrap());
}

fn search_duplications(
    strand1_file: &str,
    strand2_file: &str,

    kmer_size: usize,
    max_gap_size: u32,

    reverse: bool,
    translate: bool,
    align: bool,
    interlaced: bool,

    threads_count: usize,
    ) -> RunResult {

    let mut result : Vec<utils::SD> = Vec::new();

    let strand1 = read_fasta(strand1_file).expect(&format!("Unable to read {}", strand1_file));
    let shared_strand1 = Arc::new(strand1);

    let strand2 = {
        let mut strand2 = read_fasta(strand2_file).expect(&format!("Unable to read {}", strand2_file));
        if translate { strand2 = utils::translated(&strand2[0..strand2.len()-1].to_vec()); }
        if reverse { strand2.reverse(); }
        strand2.push(b'$');
        strand2
    };
    println!("Building suffix array...");
    let shared_suffix_array = Arc::new(r_divsufsort(&strand2));
    println!("Done.");
    let shared_strand2 = Arc::new(strand2);


    let thread_pool = ThreadPool::new(threads_count);
    let (tx, rx) = mpsc::channel();
    {
        const CHUNK_SIZE: usize = 200000;
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
                        interlaced,
                        align)).unwrap();
            });

            start += CHUNK_SIZE;
        }
    }

    drop(tx);
    println!("Looking for hulls...");
    let mut passes: Vec<utils::ProcessingSD> = rx.iter().fold(Vec::new(), |mut a, b| {a.append(&mut b.clone()); a});
    println!("Done.");


    if align {
        println!("Running perfect alignments");
        passes = passes.iter()
            .map(|b| {utils::align_perfect(b.clone())}).collect();
        println!("Done.");

        println!("Running fuzzy alignment...");
        passes = passes.iter()
            .map(|b| {utils::align_fuzzy(&(shared_strand1.clone()), &(shared_strand2.clone()), b.clone())}).collect();
        println!("Done.");
    }

    println!("Re-ordering...");
    result.extend(passes.iter().filter_map(|b| {
        match *b {
            utils::ProcessingSD::Done(ref p) => {
                Some(p.clone())
            }
            utils::ProcessingSD::ForFuzzy{ .. } => {println!("FOUND A FORFUZZY"); None}
            utils::ProcessingSD::ForSW{ .. }    => {println!("FOUND A FORSW"); None},
            utils::ProcessingSD::Empty          => {None}
        }
    }));
    println!("Done.");
    result.sort_by(|a, b|
                   if a.left != b.left {
                       (a.left).cmp(&b.left)
                   } else {
                       (a.right).cmp(&b.right)
                   });

    println!("Reducing overlapping...");
    result = reduce_overlap(&result);
    println!("Done.");

    println!("Done for {} & {}.", kmer_size, max_gap_size - kmer_size as u32);

    RunResult {
        strand1: Strand {
            length: shared_strand1.len(),
            reversed: false,
            translated: false,
        },
        strand2: Strand {
            length: shared_strand2.len() - 1, // Drop the '$'
            reversed: reverse,
            translated: translate,
        },
        SDs: result,
    }
}

// Returns true if x ⊂ y
fn subsegment((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
    let xend = xstart + xlen;
    let yend = ystart + ylen;

    xstart >= ystart && xend <= yend
}

fn overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
    let xend = xstart + xlen;
    let yend = ystart + ylen;

    (xstart >= ystart && xstart <= yend && xend >= yend)
    || (ystart >= xstart && ystart <= xend && yend >= xend)
}

fn merge(x: &utils::SD, y: &utils::SD) -> utils::SD {
   let (xleft, yleft, xsize, ysize) = if x.left < y.left {
       (x.left, y.left, x.size, y.size)
   } else {
       (y.left, x.left, x.size, y.size)
   };
   let lsize = (yleft-xleft) + ((xleft+xsize)-yleft) + ((yleft+ysize) - (xleft+xsize));

   let (xright, yright, xsize, ysize) = if x.right < y.right {
       (x.right, y.right, x.size, y.size)
   } else {
       (y.right, x.right, x.size, y.size)
   };
   let rsize = (yright-xright) + ((xright+xsize)-yright) + ((yright+ysize) - (xright+xsize));

   utils::SD {
       left: xleft,
       right: xright,
       size: cmp::min(lsize, rsize),
       rate: x.rate
   }
}

fn reduce_overlap(result: &[utils::SD]) -> Vec<utils::SD> {
    fn _reduce(result: &[utils::SD]) -> Vec<utils::SD> {
        let mut news: Vec<utils::SD> = Vec::new();
        'to_insert: for x in result.iter() {
            for ref mut y in &mut news {
                // x ⊂ y
                if subsegment(x.left_part(), y.left_part()) && subsegment(x.right_part(), y.right_part())
                {continue 'to_insert;}

                // x ⊃ y
                if subsegment(y.left_part(), x.left_part()) && subsegment(y.right_part(), x.right_part())
                {
                    y.left = x.left; y.right = x.right; y.size = x.size; y.rate = x.rate;
                    continue 'to_insert;
                }

                if overlap(x.left_part(), y.left_part()) && overlap(x.right_part(), y.right_part())
                {
                    let z = merge(x, y);
                    y.left = z.left;
                    y.right = z.right;
                    y.size = z.size;
                    continue 'to_insert;
                }
            }
            news.push(x.clone());
        }
        println!("{} reduced to {}", result.len(), news.len());
        news
    }

    let mut old_size = result.len();
    let mut news = _reduce(result);
    let mut new_size = news.len();
    while new_size < old_size {
        old_size = news.len();
        news = _reduce(&news);
        new_size = news.len();
    }
    news
}
