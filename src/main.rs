extern crate uuid;
extern crate bio;
extern crate num_cpus;
extern crate threadpool;
extern crate rustc_serialize;
extern crate env_logger;
#[macro_use]
extern crate clap;
#[macro_use]
extern crate log;
extern crate colored;
extern crate indicatif;
extern crate console;

mod automaton;
mod utils;
mod divsufsort64;
mod structs;
mod searcher;
mod logger;

use std::path;
use std::thread;
use std::time::Duration;
use std::cmp;
use std::io;
use std::io::Write;
use std::fs::File;
use std::sync::mpsc;
use std::sync::Arc;
use std::ascii::AsciiExt;
use std::time::Instant;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::TryRecvError;

use bio::io::fasta;
use log::LogLevelFilter;
use threadpool::ThreadPool;
use clap::App;
use divsufsort64::idx;
use indicatif::{ProgressBar, ProgressStyle, HumanDuration};
use console::style;

use structs::{StrandResult, RunResult, SD, Start};
use logger::Logger;

struct Strand {
    pub file_name: String,
    pub data: Vec<u8>,
    pub map: Vec<Start>,
}

fn prepare_data(strand1_file: &str,
                strand2_file: &str,
                reverse: bool,
                translate: bool,
                trim: Option<(usize, usize)>)
                -> (
                    std::sync::Arc<Strand>,
                    std::sync::Arc<Strand>,
                    std::sync::Arc<std::vec::Vec<i64>>,
                    std::sync::Arc<searcher::Searcher>,
                    ) {
    //
    // Read and map the FASTA files to process
    //
    let (map1, strand1) = read_fasta(strand1_file)
        .expect(&format!("Unable to read {}", strand1_file));
    let (map2, strand2) = {
        let (map2, strand2) = if strand2_file != strand1_file {
            read_fasta(strand2_file).expect(&format!("Unable to read {}", strand2_file))
        } else {
            trace!("Using same file for strands 1 & 2\n");
            (map1.clone(), strand1.clone())
        };
        (map2, strand2)
    };

    //
    // Ensure that shift & stop actually stay in the FASTA
    //
    let (shift, mut stop) = trim.unwrap_or((0, strand1.len() - 1));
    if stop >= strand1.len() {
        warn!("Trimming: {} greater than `{}` length ({}bp)",
              stop,
              strand1_file,
              strand1.len());
        warn!("Using {} instead of {}\n", strand1.len() - 1, stop);
        stop = strand1.len() - 1;
    }
    if stop <= shift {
        error!("ERROR: {} greater than {}", shift, stop);
        panic!();
    }
    let size = stop - shift;
    let strand1 = strand1[shift..stop].to_vec();


    //
    // Invert strands 1 & 2 to ensure efficient processing
    //
    let (mut strand1, mut strand2, _, shift2) = if strand2.len() < size {

        (
            Strand {
                file_name: strand1_file.to_owned(),
                data: strand1,
                map: map1,
            },
            Strand {
                file_name: strand2_file.to_owned(),
                data: strand2,
                map: map2,
            },
            shift,
            0
        )
    } else {
        (
            Strand {
                file_name: strand2_file.to_owned(),
                data: strand2,
                map: map2,
            },
            Strand {
                file_name: strand1_file.to_owned(),
                data: strand1,
                map: map1,
            },
            0,
            shift
        )
    };

    if translate {
        strand1.data = utils::translated(&*strand1.data);
    }
    if reverse {
        strand1.data.reverse();
    }
    strand2.data.push(b'$');

    //
    // Build the suffix array
    //
    info!("{} Building suffix array...", style("[1/4]").blue().bold());
    let suffix_array = r_divsufsort(&strand2.data);

    let shared_suffix_array = Arc::new(suffix_array);
    let shared_searcher = Arc::new(searcher::Searcher::new(&strand2.data.clone(),
                                                           &shared_suffix_array.clone(),
                                                           shift2));

    (Arc::new(strand1), Arc::new(strand2), shared_suffix_array, shared_searcher)
}

fn read_fasta(filename: &str) -> Result<(Vec<Start>, Vec<u8>), io::Error> {
    let mut map = Vec::new();
    let mut r = Vec::new();

    let reader = fasta::Reader::from_file(filename).expect("Unable to read");
    let mut counter = 0;

    for record in reader.records() {
        let record = record.expect(&format!("Unable to read {:?}: not a FASTA file", path::Path::new(filename).file_name().unwrap()));

        let name = format!("{} {}",
                           record.id(),
                           record.desc().unwrap_or(""));
        let mut seq = record.seq().to_vec();
        seq = seq.to_ascii_uppercase();
        for c in &mut seq {
            if !(structs::ALPHABET).contains(c) {
                *c = b'N'
            }
        }

        map.push(Start {
            name: name,
            position: counter,
            length: seq.len(),
        });
        counter += seq.len();
        r.append(&mut seq);
    }


    Ok((map, r))
}

pub fn r_divsufsort(dna: &[u8]) -> Vec<idx> {
    let mut sa = Vec::with_capacity(dna.len());
    sa.resize(dna.len(), 0);
    unsafe {
        divsufsort64::divsufsort64(dna.as_ptr(), sa.as_mut_ptr(), dna.len() as i64);
    }
    sa
}


fn main() {
    struct Settings {
        strand1_file: String,
        strand2_file: String,
        kmer_size: usize,
        gap_size: u32,
        min_duplication_length: usize,

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
        min_duplication_length: value_t!(args, "minlength", usize).unwrap_or_else(|_| 1000),

        reverse: args.is_present("reverse"),
        translate: args.is_present("translate"),
        interlaced: args.is_present("interlaced"),
        trim: values_t!(args, "trim", usize).unwrap_or_else(|_| Vec::new()),

        prefix: args.value_of("prefix").unwrap_or("").to_owned(),
        out: args.value_of("out").unwrap_or("").to_owned(),
        threads_count: value_t!(args, "threads", usize).unwrap_or_else(|_| num_cpus::get()),
    };

    let verbose = args.is_present("verbose");
    Logger::init(if verbose {
            LogLevelFilter::Trace
        } else {
            LogLevelFilter::Info
        })
        .unwrap();

    let out_file = if settings.out.is_empty() {
        format!("{}{}_vs_{}_{}_{}{}{}.json",
                &settings.prefix,
                path::Path::new(&settings.strand1_file).file_name().unwrap().to_str().unwrap(),
                path::Path::new(&settings.strand2_file).file_name().unwrap().to_str().unwrap(),
                settings.kmer_size,
                settings.gap_size,
                if settings.reverse {"r"} else {""},
                if settings.translate {"t"} else {""},
                )
    } else {
        settings.out + ".json"
    };

    trace!("1st strand file          {}", &settings.strand1_file);
    trace!("2nd strand file          {}", &settings.strand2_file);
    trace!("K-mers size              {}", settings.kmer_size);
    trace!("Max gap size             {}", settings.gap_size);
    trace!("Output file              {}", &out_file);
    trace!("Reverse 2nd strand       {}", settings.reverse);
    trace!("Translate 2nd strand     {}", settings.translate);
    trace!("Interlaced SD            {}", settings.interlaced);
    trace!("Threads count            {}", settings.threads_count);
    if !settings.trim.is_empty() {
        trace!("Trimming                 {} → {}\n",
              settings.trim[0],
              settings.trim[1]);
    }

    let result = search_duplications(&settings.strand1_file,
                                     &settings.strand2_file,
                                     settings.kmer_size,
                                     settings.gap_size + settings.kmer_size as u32,
                                     settings.min_duplication_length,
                                     settings.reverse,
                                     settings.translate,
                                     settings.interlaced,
                                     if !settings.trim.is_empty() {
                                         Some((settings.trim[0], settings.trim[1]))
                                     } else {
                                         None
                                     },
                                     settings.threads_count,
                                     );
    let mut out = File::create(&out_file).expect(&format!("Unable to create `{}`", &out_file));
    writeln!(&mut out,
             "{}",
             rustc_serialize::json::as_pretty_json(&result)).expect("Unable to write results");
}

#[allow(unknown_lints)]
#[allow(too_many_arguments)]
fn search_duplications(strand1_file: &str,
                       strand2_file: &str,

                       kmer_size: usize,
                       max_gap_size: u32,
                       min_duplication_length: usize,

                       reverse: bool,
                       translate: bool,
                       interlaced: bool,
                       trim: Option<(usize, usize)>,

                       threads_count: usize,
                       )
                       -> RunResult {

    let total = Instant::now();

    let (strand1, strand2, shared_suffix_array, shared_searcher) =
        prepare_data(strand1_file, strand2_file, reverse, translate, trim);


    let thread_pool = ThreadPool::new(threads_count);
    let (tx_monitor, rx_monitor) = mpsc::channel();
    let (tx, rx) = mpsc::channel();
    {
        const CHUNK_SIZE: usize = 50000;

        let num_tasks = (strand1.data.len() - kmer_size) / CHUNK_SIZE;
        let chunk_overflow = (strand1.data.len() - kmer_size) % CHUNK_SIZE;

        let mut progresses: Vec<Arc<AtomicUsize>> = Vec::with_capacity(num_tasks);

        let mut start = 0;
        for id in 0..num_tasks + 1 {
            progresses.push(Arc::new(AtomicUsize::new(0)));
            let my_tx = tx.clone();
            let suffix_array = shared_suffix_array.clone();
            let strand1 = strand1.clone();
            let strand2 = strand2.clone();
            let searcher = shared_searcher.clone();
            let my_progress = progresses[id].clone();

            thread_pool.execute(move || {
                let end = start +
                          if id < num_tasks {
                    CHUNK_SIZE
                } else {
                    chunk_overflow
                };
                my_tx.send(automaton::search_duplications(&strand1.data,
                                                         &strand2.data,
                                                         &suffix_array,
                                                         start,
                                                         end,
                                                         kmer_size,
                                                         max_gap_size,
                                                         min_duplication_length,
                                                         interlaced,
                                                         &searcher,
                                                         &my_progress))
                    .unwrap();
            });

            start += CHUNK_SIZE;
        }

        let total = strand1.data.len();
        thread::spawn(move || {
            let pb = ProgressBar::new(100);
            pb.set_style(ProgressStyle::default_bar()
                         .template("{spinner:.green} [{elapsed}] {wide_bar} {pos}% ({eta} remaining)"));

            loop {
                thread::sleep(Duration::from_millis(100));
                match rx_monitor.try_recv() {
                    Ok(_) |
                    Err(TryRecvError::Disconnected) => {
                        pb.finish_and_clear();
                        break;
                    }
                    Err(TryRecvError::Empty) => {
                        let mut current = 0;
                        for x in &progresses {
                            current += x.load(Ordering::Relaxed);
                        }
                        let percent = (current as f64 / total as f64 * 100.0) as u64;
                        pb.set_position(percent);
                    }
                }
            }
        });
    }
    drop(tx);

    info!("{} Looking for hulls...", style("[2/4]").blue().bold());
    let mut result = rx.iter().fold(Vec::new(), |mut a, b| {
        a.extend(b.iter().map(|sd| {
            SD {
                left: if !reverse {sd.left} else {strand1.data.len() - sd.left - sd.length},
                right: sd.right,
                length: sd.length,
                identity: sd.identity,
                reversed: reverse,
                translated: translate,
            }
        }));
        a
    });
    let _ = tx_monitor.send(());


    info!("{} Re-ordering...", style("[3/4]").blue().bold());
    result.sort_by(|a, b| if a.left != b.left {
        (a.left).cmp(&b.left)
    } else {
        (a.right).cmp(&b.right)
    });


    info!("{} Reducing overlapping...", style("[4/4]").blue().bold());
    result = reduce_overlap(&result);

    info!("{}",
          style(format!("{} vs. {} processed in {}.\n\n",
          strand1_file,
          strand2_file,
                  HumanDuration(total.elapsed()))).green().bold()
    );

    RunResult {
        strand1: StrandResult {
            name: strand1.file_name.to_owned(),
            length: strand1.data.len(),
            map: strand1.map.clone(),
        },
        strand2: StrandResult {
            name: strand2.file_name.to_owned(),
            length: strand2.data.len() - 1, // Drop the '$'
            map: strand2.map.clone(),
        },

        kmer: kmer_size,
        gap: max_gap_size as usize,

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

    (xstart >= ystart && xstart <= yend && xend >= yend) ||
    (ystart >= xstart && ystart <= xend && yend >= xend)
}

fn merge(x: &SD, y: &SD) -> SD {
    let (xleft, yleft, xsize, ysize) = if x.left < y.left {
        (x.left, y.left, x.length, y.length)
    } else {
        (y.left, x.left, x.length, y.length)
    };
    let lsize = (yleft - xleft) + ((xleft + xsize) - yleft) + ((yleft + ysize) - (xleft + xsize));

    let (xright, yright, xsize, ysize) = if x.right < y.right {
        (x.right, y.right, x.length, y.length)
    } else {
        (y.right, x.right, x.length, y.length)
    };
    let rsize = (yright - xright) + ((xright + xsize) - yright) +
                ((yright + ysize) - (xright + xsize));

    SD {
        left: xleft,
        right: xright,
        length: cmp::min(lsize, rsize),
        identity: x.identity,
        reversed: x.reversed,
        translated: x.translated,
    }
}

fn reduce_overlap(result: &[SD]) -> Vec<SD> {
    fn _reduce(result: &[SD]) -> Vec<SD> {
        let mut news: Vec<SD> = Vec::new();
        'to_insert: for x in result.iter() {
            for y in news.iter_mut() {
                // x ⊂ y
                if subsegment(x.left_part(), y.left_part()) &&
                   subsegment(x.right_part(), y.right_part()) {
                    continue 'to_insert;
                }

                // x ⊃ y
                if subsegment(y.left_part(), x.left_part()) &&
                   subsegment(y.right_part(), x.right_part()) {
                    y.left = x.left;
                    y.right = x.right;
                    y.length = x.length;
                    y.identity = x.identity;
                    continue 'to_insert;
                }

                if overlap(x.left_part(), y.left_part()) &&
                   overlap(x.right_part(), y.right_part()) {
                    let z = merge(x, y);
                    y.left = z.left;
                    y.right = z.right;
                    y.length = z.length;
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
