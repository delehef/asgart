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
extern crate pbr;

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
use std::time::SystemTime;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::TryRecvError;

use bio::io::fasta;
use log::LogLevelFilter;
use threadpool::ThreadPool;
use clap::App;
use divsufsort64::idx;
use pbr::ProgressBar;

use structs::{StrandResult, RunResult, SD, Start};
use logger::Logger;

struct Strand {
    pub file_name: String,
    pub data: std::sync::Arc<Vec<u8>>,
    pub map: Vec<Start>,
}

fn prepare_data(strand1_file: &str,
                strand2_file: &str,
                reverse: bool,
                translate: bool,
                trim: Option<(usize, usize)>)
                -> (Strand, Strand,
                    std::sync::Arc<std::vec::Vec<i64>>,
                    std::sync::Arc<searcher::Searcher>,
                    ) {
    //
    // Read and map the FASTA files to process
    //
    let (map1, strand1) = read_fasta(strand1_file)
        .expect(&format!("Unable to read {}", strand1_file));
    let (map2, strand2) = {
        let (map2, mut strand2) = if strand2_file != strand1_file {
            read_fasta(strand2_file).expect(&format!("Unable to read {}", strand2_file))
        } else {
            info!("Using same file for strands 1 & 2\n");
            (map1.clone(), strand1.clone())
        };
        if translate {
            strand2 = utils::translated(&strand2);
        }
        if reverse {
            strand2.reverse();
        }
        strand2.push(b'$');
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
    let (strand1, strand2, shift1, shift2) = if size > strand2.len() {
        (
            Strand {
                file_name: strand1_file.to_owned(),
                data: Arc::new(strand1),
                map: map1,
            },
            Strand {
                file_name: strand2_file.to_owned(),
                data: Arc::new(strand2),
                map: map2,
            },
            shift,
            0
        )
    } else {
        (
            Strand {
                file_name: strand2_file.to_owned(),
                data: Arc::new(strand2),
                map: map2,
            },
            Strand {
                file_name: strand1_file.to_owned(),
                data: Arc::new(strand1),
                map: map1,
            },
            0,
            shift
        )
    };

    //
    // Build the suffix array
    //
    info!("Building suffix array...");
    let now = SystemTime::now();
    let suffix_array = r_divsufsort(&strand2.data);
    info!("Done in {}s.\n", now.elapsed().unwrap().as_secs());

    let shared_suffix_array = Arc::new(suffix_array);
    let shared_searcher = Arc::new(searcher::Searcher::new(&strand2.data.clone(),
                                                           &shared_suffix_array.clone(),
                                                           shift, stop
    ));

    (strand1, strand2, shared_suffix_array, shared_searcher)
}

fn read_fasta(filename: &str) -> Result<(Vec<Start>, Vec<u8>), io::Error> {
    let mut map = Vec::new();
    let mut r = Vec::new();

    let reader = fasta::Reader::from_file(filename).expect("Unable to read");
    let mut counter = 0;

    for record in reader.records() {
        let record = record.expect(&format!("Unable to read {:?}: not a FASTA file", path::Path::new(filename).file_name().unwrap()));

        let name = format!("{} {}",
                           record.id().unwrap_or(""),
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
        min_duplication_size: usize,

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
        min_duplication_size: value_t!(args, "minsize", usize).unwrap_or_else(|_| 1000),

        reverse: args.is_present("reverse"),
        translate: args.is_present("translate"),
        interlaced: args.is_present("interlaced"),
        trim: values_t!(args, "trim", usize).unwrap_or_else(|_| Vec::new()),

        prefix: args.value_of("prefix").unwrap_or("").to_owned(),
        out: args.value_of("out").unwrap_or("").to_owned(),
        threads_count: value_t!(args, "threads", usize).unwrap_or_else(|_| num_cpus::get()),
    };

    let verbose = args.is_present("verbose");
    Logger::init(if args.is_present("verbose") {
            LogLevelFilter::Trace
        } else {
            LogLevelFilter::Off
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

    info!("1st strand file          {}", &settings.strand1_file);
    info!("2nd strand file          {}", &settings.strand2_file);
    info!("K-mers size              {}", settings.kmer_size);
    info!("Max gap size             {}", settings.gap_size);
    info!("Output file              {}", &out_file);
    info!("Reverse 2nd strand       {}", settings.reverse);
    info!("Translate 2nd strand     {}", settings.translate);
    info!("Interlaced SD            {}", settings.interlaced);
    info!("Threads count            {}", settings.threads_count);
    if !settings.trim.is_empty() {
        info!("Trimming                 {} → {}\n",
              settings.trim[0],
              settings.trim[1]);
    }

    let result = search_duplications(&settings.strand1_file,
                                     &settings.strand2_file,
                                     settings.kmer_size,
                                     settings.gap_size + settings.kmer_size as u32,
                                     settings.min_duplication_size,
                                     settings.reverse,
                                     settings.translate,
                                     settings.interlaced,
                                     if !settings.trim.is_empty() {
                                         Some((settings.trim[0], settings.trim[1]))
                                     } else {
                                         None
                                     },
                                     settings.threads_count,
                                     verbose);
    let mut out = File::create(&out_file).expect(&format!("Unable to create `{}`", &out_file));
    writeln!(&mut out,
             "{}",
             rustc_serialize::json::as_pretty_json(&result)).expect("Unable to write results");
}

fn search_duplications(strand1_file: &str,
                       strand2_file: &str,

                       kmer_size: usize,
                       max_gap_size: u32,
                       min_duplication_size: usize,

                       reverse: bool,
                       translate: bool,
                       interlaced: bool,
                       trim: Option<(usize, usize)>,

                       threads_count: usize,

                       verbose: bool)
                       -> RunResult {

    let total = SystemTime::now();

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
            let strand1 = strand1.data.clone();
            let strand2 = strand2.data.clone();
            let searcher = shared_searcher.clone();
            let my_progress = progresses[id].clone();

            thread_pool.execute(move || {
                let end = start +
                          if id < num_tasks {
                    CHUNK_SIZE
                } else {
                    chunk_overflow
                };
                my_tx.send(automaton::search_duplications(&strand1,
                                                         &strand2,
                                                         &suffix_array,
                                                         start,
                                                         end,
                                                         kmer_size,
                                                         max_gap_size,
                                                         min_duplication_size,
                                                         interlaced,
                                                         &searcher,
                                                         &my_progress))
                    .unwrap();
            });

            start += CHUNK_SIZE;
        }

        if verbose {
            let total = strand1.data.len();
            thread::spawn(move || {
                let mut pb = ProgressBar::new(100);
                loop {
                    thread::sleep(Duration::from_millis(500));
                    match rx_monitor.try_recv() {
                        Ok(_) |
                            Err(TryRecvError::Disconnected) => {
                                break;
                            }
                        Err(TryRecvError::Empty) => {
                            let mut current = 0;
                            for x in &progresses {
                                current += x.load(Ordering::Relaxed);
                            }
                            let percent = (current as f64 / total as f64 * 100.0) as u64;
                            pb.set(percent);
                        }
                    }
                }
            });
        }
    }
    drop(tx);

    info!("Looking for hulls...");
    let now = SystemTime::now();
    let mut result = rx.iter().fold(Vec::new(), |mut a, b| {
        a.extend(b.iter().map(|sd| {
            SD {
                left: sd.left,
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
    info!("Done in {}s.\n", now.elapsed().unwrap().as_secs());


    info!("Re-ordering...");
    result.sort_by(|a, b| if a.left != b.left {
        (a.left).cmp(&b.left)
    } else {
        (a.right).cmp(&b.right)
    });
    info!("Done.\n");


    info!("Reducing overlapping...");
    result = reduce_overlap(&result);
    info!("Done.\n");

    info!("{} & {} ({}/{}) processed in {}s.\n\n",
          strand1_file,
          strand2_file,
          kmer_size,
          max_gap_size - kmer_size as u32,
          total.elapsed().unwrap().as_secs());

    RunResult {
        strand1: StrandResult {
            name: path::Path::new(&strand1.file_name).file_name().unwrap().to_str().unwrap().to_owned(),
            length: strand1.data.len(),
            map: strand1.map,
        },
        strand2: StrandResult {
            name: path::Path::new(&strand2.file_name).file_name().unwrap().to_str().unwrap().to_owned(),
            length: strand2.data.len() - 1, // Drop the '$'
            map: strand2.map,
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
            for ref mut y in &mut news {
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
