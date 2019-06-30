#[macro_use] pub extern crate clap;
#[macro_use] extern crate log;
extern crate threadpool;
extern crate indicatif;
extern crate console;
extern crate num_cpus;
extern crate bio;
extern crate asgart;
extern crate rayon;
extern crate separator;

use separator::Separatable;
use indicatif::{ProgressBar, ProgressStyle, HumanDuration};
use console::style;
use bio::io::fasta;
use clap::{App, AppSettings};
use log::LevelFilter;
use rayon::prelude::*;

use std::path;
use std::thread;
use std::time::Duration;
use std::cmp;
use std::sync::mpsc;
use std::sync::Arc;
use std::time::Instant;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::TryRecvError;

use asgart::structs::*;
use asgart::logger::Logger;
use asgart::errors::*;
use asgart::exporters::Exporter;
use asgart::exporters;
use asgart::automaton;
use asgart::searcher;
use asgart::utils;
use asgart::divsufsort64;


trait Step {
    fn name(&self) -> &str;
    fn run(&self, input: Vec<ProtoSDsFamily>, strand1: &Strand, strand2: &Strand,) -> Vec<ProtoSDsFamily>;
}

struct ReOrder;
impl Step for ReOrder {
    fn name(&self) -> &str { "Re-ordering" }
    fn run(&self, mut input: Vec<ProtoSDsFamily>, strand1: &Strand, strand2: &Strand,) -> Vec<ProtoSDsFamily> {
        input
            .par_iter_mut()
            .for_each(|family| {
                family
                    .iter_mut()
                    .for_each(|sd|
                              if sd.left > sd.right && (*strand1).file_name == (*strand2).file_name {
                                  let tmp = sd.left;
                                  sd.left = sd.right;
                                  sd.right = tmp
                              }
                    )});
        input
    }
}

struct Sort;
impl Step for Sort {
    fn name(&self) -> &str { "Sorting" }
    fn run(&self, mut input: Vec<ProtoSDsFamily>, _strand1: &Strand, _strand2: &Strand,) -> Vec<ProtoSDsFamily> {
        input
            .iter_mut()
            .for_each(|family| family.sort_by(|a, b| (a.left).cmp(&b.left)));

        input
    }
}

struct ReduceOverlap;
impl Step for ReduceOverlap {
    fn name(&self) -> &str { "Reducing overlap" }
    fn run(&self, mut input: Vec<ProtoSDsFamily>, _strand1: &Strand, _strand2: &Strand,) -> Vec<ProtoSDsFamily> {
        input
            .iter_mut()
            .map(|family| reduce_overlap(&family))
            .collect()
    }
}

struct ComputeScore {trim: Option<(usize, usize)>,}
impl ComputeScore {
    fn new(trim: Option<(usize, usize)>) -> ComputeScore { ComputeScore { trim: trim.clone() }  }
}
impl Step for ComputeScore {
    fn name(&self) -> &str { "Computing Levenshtein distance" }
    fn run(&self, mut input: Vec<ProtoSDsFamily>, strand1: &Strand, strand2: &Strand) -> Vec<ProtoSDsFamily> {
        let (shift, _) = self.trim.unwrap_or((0, 0));
        input
            .par_iter_mut()
            .for_each(|family|
                      family
                      .iter_mut()
                      .for_each(|ref mut sd| sd.identity = sd.levenshtein(shift, &strand1.data, &strand2.data) as f32));
        input
    }
}

struct SearchDuplications<'a> {
    strand1: std::sync::Arc<Strand>,
    strand2: std::sync::Arc<Strand>,
    chunks_to_process: &'a[(usize, usize)],
    suffix_array: std::sync::Arc<std::vec::Vec<i64>>,
    searcher: std::sync::Arc<searcher::Searcher>,

    settings: RunSettings,
}
impl SearchDuplications<'_> {
    fn new(
        strand1: std::sync::Arc<Strand>,
        strand2: std::sync::Arc<Strand>,
        chunks_to_process: &[(usize, usize)],
        suffix_array: std::sync::Arc<std::vec::Vec<i64>>,
        searcher: std::sync::Arc<searcher::Searcher>,
        settings: RunSettings,
    ) -> SearchDuplications { SearchDuplications {
        strand1: strand1,
        strand2: strand2,
        chunks_to_process: chunks_to_process,
        suffix_array: suffix_array,
        searcher: searcher,
        settings: settings.clone(),
    }}
}
impl<'a> Step for SearchDuplications<'a> {
    fn name(&self) -> &str { "Looking for proto-duplications" }
    fn run(&self, _input: Vec<ProtoSDsFamily>, _strand1: &Strand, _strand2: &Strand) -> Vec<ProtoSDsFamily> {
        let (tx_monitor, rx_monitor) = mpsc::channel();
        let progresses = Arc::new((0..self.chunks_to_process.len()).map(|_| Arc::new(AtomicUsize::new(0))).collect::<Vec<_>>());
        let total = self.strand1.data.len();
        let monitor_thread = {
            let progresses = Arc::clone(&progresses);
            thread::spawn(move || {
                let pb = ProgressBar::new(100);
                pb.set_style(ProgressStyle::default_bar().template("{spinner:.blue} [{elapsed}] {bar:50} {pos}% ({eta} remaining)"));

                loop {
                    thread::sleep(Duration::from_millis(500));
                    match rx_monitor.try_recv() {
                        Err(TryRecvError::Empty) => {
                            pb.set_position((progresses.iter().map(|x| x.load(Ordering::Relaxed))
                                             .fold(0, |ax, x| ax + x) as f64/total as f64 * 100.0) as u64);
                        }
                        _ => { pb.finish_and_clear(); break; }
                    }
                }
            })};

        let results = self.chunks_to_process
            .par_iter()
            .enumerate()
            .map(|(id, chunk)| {
                automaton::search_duplications(
                    &self.strand1.data, &self.strand2.data,
                    &self.suffix_array, &self.searcher,
                    &progresses[id], RunSettings {start: chunk.0, end: chunk.0 + chunk.1, .. self.settings.clone()}
                )
            })
            .collect::<Vec<_>>();
        let result =
            results
            .iter()
            .fold(Vec::new(), |mut a, b| {
                a.extend(b.iter()
                         .map(|family|
                              family
                              .iter()
                              .map(|sd|{
                                  let left = if !self.settings.reverse {
                                      sd.left
                                  } else {
                                      self.strand1.data.len() - sd.left - sd.length - 1
                                  };

                                  ProtoSD {
                                      left:         left,
                                      right:        sd.right,
                                      length:       sd.length,
                                      identity:     sd.identity,
                                      reversed:     self.settings.reverse,
                                      complemented: self.settings.complement,
                                  }
                              }
                              ).collect::<ProtoSDsFamily>()));
                a
            });

        let _ = tx_monitor.send(());
        monitor_thread.join().unwrap();
        result
    }
}

type PreparedData = (
    std::sync::Arc<Strand>,                // First strand
    std::sync::Arc<Strand>,                // Second strand
    Vec<(usize, usize)>,                   // The areas that are not filled with Ns
    std::sync::Arc<std::vec::Vec<i64>>,    // Suffix array
    std::sync::Arc<searcher::Searcher>     // Searcher for the suffix array
);

struct Strand {
    pub file_name: String,
    pub data: Vec<u8>,
    pub map: Vec<Start>,
}

fn prepare_data(strand1_file: &str,
                strand2_file: &str,
                reverse: bool,
                complement: bool,
                skip_masked: bool,
                trim: Option<(usize, usize)>)
                -> Result<PreparedData> {
    fn find_chunks_to_process(strand: &[u8]) -> Vec<(usize, usize)> {
        let threshold = 5000;
        let mut start = 0;
        let mut count = 0;
        let mut r = Vec::new();

        for (i, n) in strand.iter().enumerate() {
            match n {
                b'n' | b'N' => {
                    if count != 0 {
                        if count > threshold { r.push((start, i - start)) }
                        count = 0
                    }
                }
                _ => {
                    if count == 0 {
                        count = 1;
                        start = i;
                    } else {
                        count += 1;
                    }
                }
            }
        }

        let r_length = r.iter().fold(0, |ax, c| ax + c.1);
        info!("{} chunks to process, ignoring {}bp out of {} ({}%): ",
              r.len().separated_string(),
              (strand.len() - r_length).separated_string(),
              strand.len().separated_string(),
              (((strand.len() as f64 - (r_length as f64))*100.0/strand.len() as f64) as i64).separated_string()
        );
        r.iter().for_each(|c| {
            trace!("{:>12} -> {:>12}   {:>11} bp", c.0.separated_string(), (c.0 + c.1).separated_string(), c.1.separated_string());
        });

        // We want to process the large chunks first; so that as soon as they are finished,
        // available threads can focus on processig the slower ones instead of fighting for automaton-level parallelism
        r.sort_by_key(|c| c.1);
        r
    }

    //
    // Read and map the FASTA files to process
    //
    let (map1, strand1) = read_fasta(strand1_file, skip_masked)?;
    let (map2, strand2) = {
        let (map2, strand2) = if strand2_file != strand1_file {
            read_fasta(strand2_file, skip_masked)?
        } else {
            info!("Using same file for strands 1 & 2");
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
        warn!("Using {} instead of {}", strand1.len() - 1, stop);
        stop = strand1.len() - 1;
    }
    if stop <= shift {
        return Err(format!("{} greater than {}", shift, stop).into());
    }
    let trimmed_size = stop - shift;
    let strand1 = strand1[shift..stop].to_vec();


    //
    // Invert strands 1 & 2 to ensure efficient processing
    //
    let (mut strand1, mut strand2, _, shift2) = if strand2.len() < trimmed_size {
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

    if complement { strand1.data = utils::complemented(&*strand1.data); }
    if reverse { strand1.data.reverse(); }
    strand2.data.push(b'$');

    let chunks_to_process = find_chunks_to_process(&strand1.data);

    //
    // Build the suffix array
    //
    info!("Building suffix array");
    let suffix_array = r_divsufsort(&strand2.data);

    let shared_suffix_array = Arc::new(suffix_array);
    let shared_searcher = Arc::new(
        searcher::Searcher::new(&strand2.data.clone(), &Arc::clone(&shared_suffix_array), shift2)
    );
    trace!("Done.");

    Ok((
        Arc::new(strand1), Arc::new(strand2),
        chunks_to_process,
        shared_suffix_array, shared_searcher
    ))
}

fn read_fasta(filename: &str, skip_masked: bool) -> Result<(Vec<Start>, Vec<u8>)> {
    let mut map = Vec::new();
    let mut r = Vec::new();

    let reader = fasta::Reader::from_file(filename).chain_err(|| format!("Unable to open `{}`", filename))?;
    let mut counter = 0;

    for record in reader.records() {
        let record = record.chain_err(|| format!("Unable to read {:?}: not a FASTA file", path::Path::new(filename).file_name().unwrap()))?;

        let name = record.id().to_owned();
        let mut seq = record.seq().to_vec();
        if !skip_masked {seq = seq.to_ascii_uppercase();}
        for c in &mut seq {
            if ALPHABET_MASKED.contains(c) && skip_masked {
                *c = b'N'
            } else if !(ALPHABET).contains(c) {
                trace!("Undefined base `{}` replaced by `N`", std::char::from_u32(u32::from(*c)).unwrap());
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

pub fn r_divsufsort(dna: &[u8]) -> Vec<divsufsort64::idx> {
    let mut sa = Vec::with_capacity(dna.len());
    sa.resize(dna.len(), 0);
    unsafe {
        divsufsort64::divsufsort64(dna.as_ptr(), sa.as_mut_ptr(), dna.len() as i64);
    }
    sa
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

fn merge(x: &ProtoSD, y: &ProtoSD) -> ProtoSD {
    let new_left = cmp::min(x.left, y.left);
    let lsize = cmp::max(x.left + x.length, y.left + y.length) - new_left;

    let new_right = cmp::min(x.right, y.right);
    let rsize = cmp::max(x.right + x.length, y.right + y.length) - new_right;


    ProtoSD {
        left: new_left,
        right: new_right,
        length: cmp::max(lsize, rsize),
        identity: 0.,
        reversed: x.reversed,
        complemented: x.complemented,
    }
}

fn reduce_overlap(result: &[ProtoSD]) -> Vec<ProtoSD> {
    fn _reduce(result: &[ProtoSD]) -> Vec<ProtoSD> {
        let mut news: Vec<ProtoSD> = Vec::new();
        'to_insert: for x in result.iter() {
            for y in &mut news {
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


fn run() -> Result<()> {
    // Those settings are only used to handily parse arguments
    struct Settings {
        strand1_file:           String,
        strand2_file:           String,
        kmer_size:              usize,
        gap_size:               u32,
        min_duplication_length: usize,
        max_cardinality:        usize,
        skip_masked:            bool,

        reverse:                bool,
        complement:             bool,
        trim:                   Vec<usize>,

        prefix:                 String,
        out:                    String,
        out_format:             String,
        compute_score:          bool,
        threads_count:          usize,
    }

    let yaml = load_yaml!("asgart.yaml");
    let args = App::from_yaml(yaml)
        .version(crate_version!())
        .author(crate_authors!())
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .get_matches();

    let settings = Settings {
        strand1_file:           args.value_of("strand1").unwrap().to_owned(),
        strand2_file:           args.value_of("strand2").unwrap().to_owned(),
        kmer_size:              value_t_or_exit!(args, "probe_size", usize),
        gap_size:               value_t_or_exit!(args, "max_gap", u32),
        min_duplication_length: value_t!(args, "min_length", usize).unwrap(),
        max_cardinality:        value_t!(args, "max_cardinality", usize).unwrap(),
        skip_masked:            args.is_present("skip_masked"),

        reverse:                args.is_present("reverse"),
        complement:             args.is_present("complement"),
        trim:                   values_t!(args, "trim", usize).unwrap_or_else(|_| Vec::new()),

        prefix:                 args.value_of("prefix").unwrap().to_owned(),
        out:                    args.value_of("out").unwrap_or("").to_owned(),
        out_format:             args.value_of("out_format").unwrap().to_owned(),
        compute_score:          args.is_present("compute_score"),
        threads_count:          value_t!(args, "threads", usize).unwrap_or_else(|_| num_cpus::get()),
    };

    Logger::init(if args.is_present("verbose") {LevelFilter::Trace} else {LevelFilter::Info}).unwrap();

    let out_file = if settings.out.is_empty() {
        format!("{}{}-{}{}{}{}",
                &settings.prefix,
                path::Path::new(&settings.strand1_file).file_stem().unwrap().to_str().unwrap(),
                path::Path::new(&settings.strand2_file).file_stem().unwrap().to_str().unwrap(),
                if settings.reverse || settings.complement {"_"} else {""},
                if settings.reverse {"R"} else {""},
                if settings.complement {"C"} else {""},
        )
    } else {
        settings.out
    };

    info!("1st strand file          {}", &settings.strand1_file);
    info!("2nd strand file          {}", &settings.strand2_file);
    info!("K-mers size              {}", settings.kmer_size);
    info!("Max gap size             {}", settings.gap_size);
    info!("Output file              {}", &out_file);
    info!("Reverse 2nd strand       {}", settings.reverse);
    info!("Complement 2nd strand    {}", settings.complement);
    info!("Skipping soft-masked     {}", settings.skip_masked);
    info!("Min. length              {}", settings.min_duplication_length);
    info!("Max. cardinality         {}", settings.max_cardinality);
    info!("Threads count            {}", settings.threads_count);
    if !settings.trim.is_empty() {
        info!("Trimming                   {} → {}", settings.trim[0], settings.trim[1]);
    }


    let result = search_duplications(
        &settings.strand1_file,
        &settings.strand2_file,
        RunSettings {
            probe_size:             settings.kmer_size,
            max_gap_size:           settings.gap_size + settings.kmer_size as u32,
            min_duplication_length: settings.min_duplication_length,
            max_cardinality:        settings.max_cardinality,

            reverse:                settings.reverse,
            complement:             settings.complement,
            skip_masked:            settings.skip_masked,

            start:                  0,
            end:                    0,

            compute_score:          settings.compute_score,
            threads_count:          settings.threads_count,
            trim:                   if !settings.trim.is_empty() {
                Some((settings.trim[0], settings.trim[1]))
            } else {
                None
            },
        },
    )?;

    let exporter = match &settings.out_format[..] {
        "json" => { Box::new(exporters::JSONExporter) as Box<dyn Exporter> }
        "gff2" => { Box::new(exporters::GFF2Exporter) as Box<dyn Exporter> }
        "gff3" => { Box::new(exporters::GFF3Exporter) as Box<dyn Exporter> }
        _      => { Box::new(exporters::JSONExporter) as Box<dyn Exporter> }
    };
    let out_file_name = exporter.save(&result, &out_file)?;
    info!("{}", style(format!("Result written to {}", &out_file_name)).bold());
    Ok(())
}


fn search_duplications(
    strand1_file: &str,
    strand2_file: &str,
    settings: RunSettings,
) -> Result<RunResult> {

    let total = Instant::now();

    info!("Preprocessing data");
    let (strand1, strand2, to_process, shared_suffix_array, shared_searcher) =
        prepare_data(
            strand1_file,
            strand2_file,
            settings.reverse,
            settings.complement,
            settings.skip_masked,
            settings.trim)?;



    let mut steps : Vec<Box<dyn Step>> = Vec::new();
    steps.push(Box::new(SearchDuplications::new(
        strand1.clone(),
        strand2.clone(),
        &to_process,
        shared_suffix_array,
        shared_searcher,
        settings,
    )));
    if settings.compute_score {steps.push(Box::new(ComputeScore::new(settings.trim)));}
    steps.push(Box::new(ReOrder{}));
    steps.push(Box::new(ReduceOverlap{}));

    let mut result = Vec::new();
    for (i, step) in steps.iter().enumerate() {
        info!("{} {}...", style(format!("[{}/{}]", i + 1, steps.len())).blue().bold(), step.name());
        result = step.run(result, &strand1, &strand2);
    }

    info!("{}",
          style(format!("{:?} vs. {:?} processed in {}.",
                        path::Path::new(strand1_file).file_name().expect("Should never fail"),
                        path::Path::new(strand2_file).file_name().expect("Should never fail"),
                        HumanDuration(total.elapsed()))).green().bold()
    );

    let strand1 = StrandResult {
        name: strand1.file_name.to_owned(),
        length: strand1.map.iter().fold(0, |ax, chr| ax + chr.length),
        map: strand1.map.clone(),
    };
    let strand2 = StrandResult {
        name: strand2.file_name.to_owned(),
        length: strand2.map.iter().fold(0, |ax, chr| ax + chr.length),
        map: strand2.map.clone(),
    };

    Ok(RunResult {
        strand1: strand1.clone(),
        strand2: strand2.clone(),
        settings: settings,
        families: result.iter()
            .map(|family|
                 family.iter()
                 .map(|sd|
                      {
                          SD {
                              chr_left: strand1.find_chr_by_pos(sd.left).name.clone(),
                              chr_right: strand2.find_chr_by_pos(sd.right).name.clone(),

                              global_left_position: sd.left,
                              global_right_position: sd.right,

                              chr_left_position: sd.left - strand1.find_chr_by_pos(sd.left).position,
                              chr_right_position: sd.right - strand2.find_chr_by_pos(sd.right).position,

                              length: sd.length,
                              identity: sd.identity,
                              reversed: sd.reversed,
                              complemented: sd.complemented,
                          }})
                 .collect::<Vec<SD>>()
            ).collect(),
    })
}


fn main() {
    if let Err(ref e) = run() {
        println!("{} {}", style("Error: ").red(), e);
        for e in e.iter().skip(1) {
            println!("{}", e);
        }
        if let Some(backtrace) = e.backtrace() {
            println!("backtrace: {:?}", backtrace);
        }
        std::process::exit(1);
    }
}
