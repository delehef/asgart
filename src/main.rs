extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;
extern crate rustc_serialize;
#[macro_use]
extern crate clap;

use std::{thread, time};
use std::time::Duration;
use std::cmp;
use std::io;
use std::io::BufRead;
use std::io::Write;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;
use std::ascii::AsciiExt;
use std::time::SystemTime;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::TryRecvError;

use threadpool::ThreadPool;

use clap::App;

use divsufsort64::idx;

mod utils;
mod divsufsort64;
mod structs;
mod searcher;

use structs::{Strand,RunResult,SD};

static mut VERBOSE: bool = false;

fn read_fasta(filename: &str) -> Result<Vec<u8>, io::Error> {
    let file = try!(File::open(filename));
    let file = io::BufReader::new(file);


    let mut r = file.lines()
       .filter_map(|result| result.ok())
       .filter(|line| !line.starts_with('>'))
       .fold(Vec::new(), |mut r, line| {r.extend(line.trim().as_bytes().iter().cloned()); r});

    r = r.to_ascii_uppercase();

    // r.retain(|c| *c != b'n' && *c != b'N');
    r.retain(|c| *c != b'>');

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


macro_rules! log(
    ($($arg:tt)*) => (
        let v = unsafe { VERBOSE };
        if v {
            match writeln!(&mut ::std::io::stderr(), $($arg)* ) {
                Ok(_) => {},
                Err(_) => {}
            }
        }
        )
);

fn main() {
    struct Settings {
        strand1_file: String,
        strand2_file: String,
        kmer_size: usize,
        gap_size: u32,

        reverse: bool,
        translate: bool,
        interlaced: bool,
        trim: Vec<usize>,

        prefix: String,
        out: String,
        threads_count: usize,
    }

    let yaml = load_yaml!("cli.yaml");
    let args = App::from_yaml(yaml)
        .version(crate_version!())
        .author(crate_authors!())
        .get_matches();

    let settings = Settings {
        strand1_file: args.value_of("strand1").unwrap().to_owned(),
        strand2_file: args.value_of("strand2").unwrap().to_owned(),
        kmer_size: value_t_or_exit!(args, "probe_size", usize),
        gap_size: value_t_or_exit!(args, "max_gap", u32),

        reverse: args.is_present("reverse"),
        translate: args.is_present("translate"),
        interlaced: args.is_present("interlaced"),
        trim: values_t!(args.values_of("trim"), usize).unwrap_or_else(|_| Vec::new()),

        prefix: args.value_of("prefix").unwrap_or("").to_owned(),
        out: args.value_of("out").unwrap_or("").to_owned(),
        threads_count: value_t!(args, "threads", usize).unwrap_or_else(|_| num_cpus::get()),
        
    };

    unsafe { VERBOSE = args.is_present("verbose"); }

    let out_file = if settings.out.is_empty() {
        format!("{}{}_vs_{}_{}_{}{}{}.json",
                &settings.prefix,
                &settings.strand1_file,
                &settings.strand2_file,
                settings.kmer_size,
                settings.gap_size,
                if settings.reverse {"r"} else {""},
                if settings.translate {"t"} else {""},
                )
    } else {
        settings.out + ".json"
    };

    if unsafe {VERBOSE} {
        println!("1st strand file          {}", &settings.strand1_file);
        println!("2nd strand file          {}", &settings.strand2_file);
        println!("K-mers size              {}", settings.kmer_size);
        println!("Max gap size             {}", settings.gap_size);
        println!("Output file              {}", &out_file);
        println!("Reverse 2nd strand       {}", settings.reverse);
        println!("Translate 2nd strand     {}", settings.translate);
        println!("Interlaced SD            {}", settings.interlaced);
        println!("Threads count            {}", settings.threads_count);
        if settings.trim.len() > 0 {
            println!("Trimming                 {} - {}", settings.trim[0], settings.trim[1]);
        }
        println!("libdivsufsort            v{:?}", unsafe{divsufsort64::divsufsort64_version()});
        println!("");
    }

    let result = search_duplications(
        &settings.strand1_file, &settings.strand2_file,
        settings.kmer_size, settings.gap_size + settings.kmer_size as u32,
        settings.reverse, settings.translate, settings.interlaced,
        if !settings.trim.is_empty() {Some((settings.trim[0], settings.trim[1]))} else {None},
        settings.threads_count,
        );
    let mut out = File::create(&out_file).expect(&format!("Unable to create `{}`", &out_file));
    writeln!(&mut out, "{}", rustc_serialize::json::encode(&result).unwrap());
}

fn search_duplications(
    strand1_file: &str,
    strand2_file: &str,

    kmer_size: usize,
    max_gap_size: u32,

    reverse: bool,
    translate: bool,
    interlaced: bool,
    trim: Option<(usize, usize)>,

    threads_count: usize,) -> RunResult {
    let total = SystemTime::now();

    let strand1 = read_fasta(strand1_file).expect(&format!("Unable to read {}", strand1_file));
    let shared_strand1 = Arc::new(strand1);

    let strand2 = {
        let mut strand2 = read_fasta(strand2_file).expect(&format!("Unable to read {}", strand2_file));
        if translate { strand2 = utils::translated(&strand2[0..strand2.len()-1].to_vec()); }
        if reverse { strand2.reverse(); }
        strand2.push(b'$');
        strand2
    };

    let (shift, stop) = trim.unwrap_or((0, shared_strand1.len()));
    if stop > shared_strand1.len() {
        panic!("ERROR: {} greater than `{}` length ({}bp)", stop, strand1_file,
               shared_strand1.len());
    }
    if stop <= shift {
        panic!("ERROR: {} greater than {}", shift, stop);
    }

    log!("Building suffix array...");
    let now = SystemTime::now();
    let shared_suffix_array = Arc::new(r_divsufsort(&strand2));
    let shared_searcher = Arc::new(searcher::Searcher::new(&strand2, 
                                                           &shared_suffix_array.clone()));
    let shared_strand2 = Arc::new(strand2);
    log!("Done in {}s.", now.elapsed().unwrap().as_secs());


    let thread_pool = ThreadPool::new(threads_count);
    let (tx_monitor, rx_monitor) = mpsc::channel();
    let (tx, rx) = mpsc::channel();
    {
        const CHUNK_SIZE: usize = 200000;
        let size = stop - shift;

        let num_tasks = (size - kmer_size)/CHUNK_SIZE;
        let chunk_overflow = (size - kmer_size)%CHUNK_SIZE;

        let mut progresses: Vec<Arc<AtomicUsize>> = Vec::with_capacity(num_tasks);

        let mut start = shift;
        for id in 0..num_tasks+1
        {
            progresses.push(Arc::new(AtomicUsize::new(0)));
            let suffix_array = shared_suffix_array.clone();
            let strand1 = shared_strand1.clone();
            let strand2 = shared_strand2.clone();
            let searcher = shared_searcher.clone();

            let my_tx = tx.clone();
            let my_progress = progresses[id].clone();

            thread_pool.execute(move || {
                let end = start + if id<num_tasks {CHUNK_SIZE} else {chunk_overflow};
                my_tx.send(utils::search_duplications(
                        &strand1, &strand2, &suffix_array,
                        start, end,
                        kmer_size, max_gap_size,
                        interlaced, &searcher,
                        &my_progress)).unwrap();
            });

            start += CHUNK_SIZE;
        }

        let total = shared_strand1.len();
        thread::spawn(move || {
            loop {
                thread::sleep(Duration::from_millis(1000));
                match rx_monitor.try_recv() {
                    Ok(_) | Err(TryRecvError::Disconnected) => {
                        break;
                    }
                    Err(TryRecvError::Empty) => {
                        let mut current = 0;
                        for x in progresses.iter() {
                            current += x.load(Ordering::Relaxed);
                        }
                        log!("{}%", (current as f64/total as f64 * 100.0) as u32);
                    }
                }
            }
        });
    }
    drop(tx);

    log!("Looking for hulls...");
    let now = SystemTime::now();
    let mut result = rx.iter().fold(Vec::new(), |mut a, b| {a.append(&mut b.clone()); a});
    let _ = tx_monitor.send(());
    log!("Done in {}s.", now.elapsed().unwrap().as_secs());


    log!("Re-ordering...");
    result.sort_by(|a, b|
                   if a.left != b.left {
                       (a.left).cmp(&b.left)
                   } else {
                       (a.right).cmp(&b.right)
                   });
    log!("Done.");


    log!("Reducing overlapping...");
    result = reduce_overlap(&result);
    log!("Done.");

    log!("{} & {} ({}/{}) processed in {}s.",
        strand1_file, strand2_file, kmer_size, max_gap_size - kmer_size as u32,
        total.elapsed().unwrap().as_secs());

    RunResult {
        strand1: Strand {
            name: strand1_file.to_owned(),
            length: shared_strand1.len(),
            reversed: false,
            translated: false,
        },
        strand2: Strand {
            name: strand2_file.to_owned(),
            length: shared_strand2.len() - 1, // Drop the '$'
            reversed: reverse,
            translated: translate,
        },
        sds: result,
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

fn merge(x: &SD, y: &SD) -> SD {
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

   SD {
       left: xleft,
       right: xright,
       size: cmp::min(lsize, rsize),
       identity: x.identity
   }
}

fn reduce_overlap(result: &[SD]) -> Vec<SD> {
    fn _reduce(result: &[SD]) -> Vec<SD> {
        let mut news: Vec<SD> = Vec::new();
        'to_insert: for x in result.iter() {
            for ref mut y in &mut news {
                // x ⊂ y
                if subsegment(x.left_part(), y.left_part()) && subsegment(x.right_part(), y.right_part())
                {continue 'to_insert;}

                // x ⊃ y
                if subsegment(y.left_part(), x.left_part()) && subsegment(y.right_part(), x.right_part())
                {
                    y.left = x.left; y.right = x.right; y.size = x.size; y.identity = x.identity;
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
