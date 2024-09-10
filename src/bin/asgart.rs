use std::{
    cmp, path,
    sync::{
        atomic::{AtomicUsize, Ordering},
        mpsc::{self, TryRecvError},
        Arc,
    },
    thread,
    time::{Duration, Instant},
};

use anyhow::{Context, Result};
use clap::*;
use console::style;
use indicatif::{HumanDuration, ProgressBar, ProgressStyle};
use log::*;
use rayon::prelude::*;
use thousands::Separable;

use asgart::{
    automaton,
    divsufsort::{divsufsort64, SuffixArray},
    exporters,
    logger::Logger,
    searcher,
    structs::*,
    utils,
};

trait Step {
    fn name(&self) -> &str;
    fn run(&self, input: Vec<ProtoSDsFamily>, strand: &Strand) -> Vec<ProtoSDsFamily>;
}

struct ReOrder;
impl Step for ReOrder {
    fn name(&self) -> &str {
        "Re-ordering"
    }

    fn run(&self, mut input: Vec<ProtoSDsFamily>, _strand: &Strand) -> Vec<ProtoSDsFamily> {
        input.par_iter_mut().for_each(|family| {
            family.iter_mut().for_each(|sd| {
                if sd.left > sd.right {
                    let tmp = sd.left;
                    sd.left = sd.right;
                    sd.right = tmp
                }
            })
        });
        input
    }
}

struct Sort;
impl Step for Sort {
    fn name(&self) -> &str {
        "Sorting"
    }

    fn run(&self, mut input: Vec<ProtoSDsFamily>, _strand: &Strand) -> Vec<ProtoSDsFamily> {
        input
            .iter_mut()
            .for_each(|family| family.sort_by(|a, b| (a.left).cmp(&b.left)));
        input
    }
}

struct ReduceOverlap;
impl Step for ReduceOverlap {
    fn name(&self) -> &str {
        "Reducing overlap"
    }

    fn run(&self, mut input: Vec<ProtoSDsFamily>, _strand: &Strand) -> Vec<ProtoSDsFamily> {
        input
            .iter_mut()
            .map(|family| reduce_overlap(&family))
            .collect()
    }
}

struct FilterNs;
impl Step for FilterNs {
    fn name(&self) -> &str {
        "Filtering uncertain duplications"
    }

    fn run(&self, mut input: Vec<ProtoSDsFamily>, strand: &Strand) -> Vec<ProtoSDsFamily> {
        input
            .par_iter_mut()
            .for_each(|family| family.retain(|sd| sd.n_content(&strand.data) <= 0.2));
        input
            .into_iter()
            .filter(|family| !family.is_empty())
            .collect::<Vec<_>>()
    }
}

struct ComputeScore;
impl Step for ComputeScore {
    fn name(&self) -> &str {
        "Computing Levenshtein distance"
    }

    fn run(&self, mut input: Vec<ProtoSDsFamily>, strand: &Strand) -> Vec<ProtoSDsFamily> {
        input.par_iter_mut().for_each(|family| {
            family
                .iter_mut()
                .for_each(|ref mut sd| sd.identity = sd.levenshtein(&strand.data) as f32)
        });
        input
    }
}

struct SearchDuplications<'a> {
    chunks_to_process: &'a [(usize, usize)],
    trim: Option<(usize, usize)>,
    settings: RunSettings,
}
impl SearchDuplications<'_> {
    fn new(
        chunks_to_process: &[(usize, usize)],
        trim: Option<(usize, usize)>,
        settings: RunSettings,
    ) -> SearchDuplications {
        SearchDuplications {
            chunks_to_process: chunks_to_process,
            trim: trim,
            settings: settings.clone(),
        }
    }
}
impl<'a> Step for SearchDuplications<'a> {
    fn name(&self) -> &str {
        "Looking for proto-duplications"
    }

    fn run(&self, _input: Vec<ProtoSDsFamily>, strand: &Strand) -> Vec<ProtoSDsFamily> {
        // Build the suffix array
        //
        debug!("Building suffix array");
        let sa_build_time = Instant::now();
        let shared_suffix_array = if let Some((start, end)) = self.trim {
            let mut sub_strand = strand.data[start..end].to_vec();
            sub_strand.push(b'$');
            let mut suffix_array = r_divsufsort(&sub_strand);
            suffix_array.iter_mut().for_each(|x| *x += start as i64);
            Arc::new(suffix_array)
        } else {
            Arc::new(r_divsufsort(&strand.data))
        };
        let shared_searcher = Arc::new(searcher::Searcher::new(
            &strand.data,
            &Arc::clone(&shared_suffix_array),
            0,
        ));
        debug!("Done in {}", HumanDuration(sa_build_time.elapsed()));

        // Set up th progress bar
        //
        let (tx_monitor, rx_monitor) = mpsc::channel();
        let progresses = Arc::new(
            (0..self.chunks_to_process.len())
                .map(|_| Arc::new(AtomicUsize::new(0)))
                .collect::<Vec<_>>(),
        );
        let total = self.chunks_to_process.iter().fold(0, |ax, c| ax + c.1);
        let monitor_thread = {
            let progresses = Arc::clone(&progresses);
            thread::spawn(move || {
                let pb = ProgressBar::new(100);
                pb.set_style(
                    ProgressStyle::default_bar()
                        .template("{spinner:.blue} [{elapsed}] {bar:50} {pos}% (~{eta} remaining)")
                        .unwrap(),
                );

                loop {
                    thread::sleep(Duration::from_millis(500));
                    match rx_monitor.try_recv() {
                        Err(TryRecvError::Empty) => {
                            pb.set_position(
                                (progresses
                                    .iter()
                                    .map(|x| x.load(Ordering::Relaxed))
                                    .fold(0, |ax, x| ax + x)
                                    as f64
                                    / total as f64
                                    * 100.0) as u64,
                            );
                        }
                        _ => {
                            pb.finish_and_clear();
                            break;
                        }
                    }
                }
            })
        };

        // And do the job
        //
        let results = self
            .chunks_to_process
            .par_iter()
            .enumerate()
            .map(|(id, chunk)| {
                let mut _needle;
                let needle = if !self.settings.reverse && !self.settings.complement {
                    &strand.data[chunk.0..chunk.0 + chunk.1]
                } else {
                    _needle = strand.data[chunk.0..chunk.0 + chunk.1].to_vec();
                    if self.settings.complement {
                        _needle = utils::complemented(&_needle);
                    }
                    if self.settings.reverse {
                        _needle.reverse();
                    }
                    &_needle
                };

                let mut proto_sds_families = automaton::search_duplications(
                    id,
                    needle,
                    chunk.0,
                    &strand.data,
                    &shared_suffix_array.clone(),
                    &shared_searcher.clone(),
                    &progresses[id],
                    self.settings.clone(),
                );
                proto_sds_families.iter_mut().for_each(|proto_family| {
                    proto_family.iter_mut().for_each(|proto_sd| {
                        if !self.settings.reverse {
                            proto_sd.left += chunk.0
                        } else {
                            proto_sd.left = chunk.0 + chunk.1 - proto_sd.left - proto_sd.left_length
                        }
                    })
                });
                proto_sds_families
            })
            .collect::<Vec<_>>();
        let result = results.iter().fold(Vec::new(), |mut a, b| {
            a.extend(b.iter().map(|family| {
                family
                    .into_iter()
                    .map(|sd| ProtoSD {
                        reversed: self.settings.reverse,
                        complemented: self.settings.complement,
                        ..*sd
                    })
                    .collect::<ProtoSDsFamily>()
            }));
            a
        });

        let _ = tx_monitor.send(());
        monitor_thread.join().unwrap();
        result
    }
}

type PreparedData = (
    Option<(usize, usize)>, // Trim
    Vec<(usize, usize)>,    // The areas that are not filled with Ns
    Strand,                 // DNA strand to process
);

struct Strand {
    pub file_names: String,
    pub data: Vec<u8>,
    pub map: Vec<Start>,
}

fn prepare_data(
    strands_files: &[String],
    skip_masked: bool,
    trim: Option<(usize, usize)>,
) -> Result<PreparedData> {
    fn read_fasta(filename: &str, skip_masked: bool) -> Result<(Vec<Start>, Vec<u8>)> {
        let mut map = Vec::new();
        let mut r = Vec::new();

        let reader = bio::io::fasta::Reader::from_file(filename)
            .with_context(|| format!("Unable to read FASTA file `{}`", filename))?;
        let mut counter = 0;

        for record in reader.records() {
            let record = record.context(format!("Unable to parse `{}`", filename))?;

            let name = record.id().to_owned();
            let mut seq = record.seq().to_vec();
            if !skip_masked {
                seq = seq.to_ascii_uppercase();
            }
            for c in &mut seq {
                if ALPHABET_MASKED.contains(c) && skip_masked {
                    *c = b'N'
                } else if !(ALPHABET).contains(c) {
                    trace!("Undefined base `{}` replaced by `N`", *c as char);
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

    // Given a DNA fragment, returns a list of segments without too many N's in them
    // Coordinates are relative to the fragment.
    fn find_chunks_to_process(strand: &[u8]) -> Vec<(usize, usize)> {
        fn count_n(strand: &[u8], start: usize) -> usize {
            strand
                .iter()
                .skip(start)
                .take_while(|x| **x == b'n' || **x == b'N')
                .count()
        }

        let threshold = 5000;
        let mut start = 0;
        let mut count = 0;
        let mut chunks = Vec::new();
        let mut i = 0;
        while i < strand.len() {
            let n = strand[i];
            match n {
                b'n' | b'N' => {
                    let n_count = count_n(&strand, i);
                    if n_count > threshold {
                        if count > 0 {
                            chunks.push((start, count));
                            count = 0;
                        };
                        start = i + n_count;
                    } else {
                        count += n_count;
                    }
                    i += n_count;
                }
                _ => {
                    if count == 0 {
                        count = 1;
                        start = i;
                    } else {
                        count += 1;
                    }
                    i += 1
                }
            }
        }
        if count != 0 {
            chunks.push((start, count))
        };
        if chunks.is_empty() {
            chunks.push((0, strand.len()))
        };

        chunks
    }

    // Read and map the FASTA files to process
    //
    let mut maps = Vec::new();
    let mut strand = Vec::new();
    let mut offset = 0;
    let mut chunks_to_process = Vec::new();

    for file_name in strands_files {
        let (map, new_strand) = read_fasta(file_name, skip_masked)
            .with_context(|| format!("Unable to parse `{}`", file_name))?;

        // We want to add each fragment separately to ensure that chunks are cutting
        // between fragments
        for chr in map.iter() {
            chunks_to_process.extend(
                find_chunks_to_process(&new_strand[chr.position..chr.position + chr.length])
                    .into_iter()
                    .map(|(start, length)| (chr.position + offset + start, length)),
            );
        }
        maps.extend(map.into_iter().map(|start| Start {
            position: start.position + offset,
            ..start
        }));

        offset = offset + new_strand.len();
        strand.extend(new_strand);
    }
    info!(
        "Parsed {} file{} containing a total of {} fragments",
        strands_files.len(),
        if strands_files.len() > 1 { "s" } else { "" },
        maps.len()
    );
    maps.iter().for_each(|s| {
        debug!(
            "{:>20}: {:>15}  --> {:>15}    {:>15} bp",
            s.name,
            s.position.separate_with_spaces(),
            (s.position + s.length).separate_with_spaces(),
            s.length.separate_with_spaces()
        )
    });

    let chunks_length = chunks_to_process.iter().fold(0, |ax, c| ax + c.1);
    info!(
        "Processing {} chunks totalling {}bp, skipping {}bp out of {} ({}%)",
        chunks_to_process.len().separate_with_spaces(),
        chunks_length.separate_with_spaces(),
        (strand.len() - chunks_length).separate_with_spaces(),
        strand.len().separate_with_spaces(),
        (((strand.len() as f64 - (chunks_length as f64)) * 100.0 / strand.len() as f64) as i64)
            .separate_with_spaces()
    );
    chunks_to_process.iter().for_each(|c| {
        trace!(
            "{:>12} -> {:>12}   {:>11} bp",
            c.0.separate_with_spaces(),
            (c.0 + c.1).separate_with_spaces(),
            c.1.separate_with_spaces()
        );
    });
    strand.push(b'$'); // For the SA construction

    Ok((
        trim.and_then(|(shift, _stop)| {
            // Ensure that shift & stop actually stay in the dataset
            //
            let mut stop = _stop;
            if stop >= strand.len() {
                warn!(
                    "Trimming: {} greater than total length ({}bp)",
                    stop,
                    strand.len()
                );
                warn!("Using {} instead of {}", strand.len() - 1, stop);
                stop = strand.len() - 1;
            }

            if stop <= shift {
                warn!(
                    "Trimming: {} greater than {}, skipping trimming",
                    shift, stop
                );
                None
            } else if shift >= strand.len() {
                warn!(
                    "Trimming: {} greater than total length ({}bp), skipping trimming",
                    shift,
                    strand.len()
                );
                None
            } else {
                Some((shift, stop))
            }
        }),
        chunks_to_process,
        Strand {
            file_names: strands_files.join(", "),
            data: strand,
            map: maps,
        },
    ))
}

pub fn r_divsufsort(dna: &[u8]) -> SuffixArray {
    let mut sa = Vec::with_capacity(dna.len());
    sa.resize(dna.len(), 0);
    unsafe {
        divsufsort64(dna.as_ptr(), sa.as_mut_ptr(), dna.len() as i64);
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

    (xstart >= ystart && xstart <= yend && xend >= yend)
        || (ystart >= xstart && ystart <= xend && yend >= xend)
}

fn merge(x: &ProtoSD, y: &ProtoSD) -> ProtoSD {
    let new_left = cmp::min(x.left, y.left);
    let lsize = cmp::max(x.left + x.left_length, y.left + y.right_length) - new_left;

    let new_right = cmp::min(x.right, y.right);
    let rsize = cmp::max(x.right + x.left_length, y.right + y.right_length) - new_right;

    ProtoSD {
        left: new_left,
        right: new_right,
        left_length: lsize,
        right_length: rsize,
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
                if subsegment(x.left_part(), y.left_part())
                    && subsegment(x.right_part(), y.right_part())
                {
                    continue 'to_insert;
                }

                // x ⊃ y
                if subsegment(y.left_part(), x.left_part())
                    && subsegment(y.right_part(), x.right_part())
                {
                    y.left = x.left;
                    y.right = x.right;
                    y.left_length = x.left_length;
                    y.right_length = x.right_length;
                    continue 'to_insert;
                }

                if overlap(x.left_part(), y.left_part()) && overlap(x.right_part(), y.right_part())
                {
                    let z = merge(x, y);
                    y.left = z.left;
                    y.right = z.right;
                    y.left_length = z.left_length;
                    y.right_length = z.right_length;
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

#[derive(Parser)]
#[command(
    name = "ASGART",
    version,
    author,
    about = "A Segmental duplications Gathering and Refinement Tool"
)]
struct Args {
    #[arg()]
    strands: Vec<String>,

    #[arg(
        long,
        help = "minimal length (in bp) of the duplications",
        default_value = "1000"
    )]
    min_length: usize,

    #[arg(
        short = 'k',
        long,
        help = "length of the probing k-mers",
        default_value = "20"
    )]
    probe_size: usize,

    #[arg(short = 'g', long, default_value = "100")]
    /// Maximum length of a gap
    gap_size: usize,

    #[arg(short = 'R', long)]
    /// Search for reversed duplications
    reverse: bool,

    #[arg(short = 'C', long)]
    /// Search for complemented duplications
    complement: bool,

    #[arg(short = 'S', long)]
    /// Ignore soft-masked repeated zones (lowercased regions)
    skip_masked: bool,

    #[arg(long, num_args = 2)]
    /// Trim the first strand
    trim: Option<Vec<usize>>,

    #[arg(long, default_value = "500")]
    /// maximal cardinality of duplication families
    max_cardinality: usize,

    #[arg(long, default_value = "")]
    /// prefix to prepend to the default output file name
    prefix: String,

    #[arg(long)]
    /// set the output file name
    out: Option<String>,

    #[arg(long)]
    /// Compute the Levenshtein distance between duplicons
    /// /!\ WARNING THIS IS A TIME- AND MEMORY-HEAVY OPERATION
    compute_score: bool,

    #[arg(long)]
    /// number of threads to use; default to the number of cores
    threads: Option<usize>,

    #[arg(long, default_value = "1000000")]
    /// Size used to slice input data for parallel processing
    chunk_size: usize,
}

fn main() -> Result<()> {
    let args = Args::parse();

    // Logger::init(match args.occurrences_of("verbose") {
    //     0 => LevelFilter::Info,
    //     1 => LevelFilter::Debug,
    //     2 => LevelFilter::Trace,
    //     _ => LevelFilter::Trace,
    // })
    // .with_context(|| "Unable to initialize logger")?;

    let radix = (args
        .strands
        .iter()
        .map(|n| {
            path::Path::new(&n)
                .file_stem()
                .expect(&format!("Unable to parse {}", n))
                .to_str()
                .unwrap()
                .to_string()
        })
        .collect::<Vec<String>>())
    .join("-");

    info!("Processing {}", &args.strands.join(", "));
    debug!("K-mers size                {}", args.probe_size);
    debug!("Max gap size               {}", args.gap_size);
    debug!("Reversed duplications      {}", args.reverse);
    debug!("Complemented duplications  {}", args.complement);
    debug!("Skipping soft-masked       {}", args.skip_masked);
    debug!("Min. length                {}", args.min_length);
    debug!("Max. cardinality           {}", args.max_cardinality);
    debug!(
        "Threads count              {}",
        args.threads.unwrap_or(num_cpus::get_physical())
    );
    if let Some(trim) = args.trim.as_ref() {
        debug!("Trimming                   {} → {}", trim[0], trim[1]);
    }

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads.unwrap_or(num_cpus::get_physical()))
        .build_global()
        .expect("Unable to create thread pool");

    let result = search_duplications(
        &args.strands,
        RunSettings {
            probe_size: args.probe_size,
            max_gap_size: args.gap_size as u32 + args.probe_size as u32,
            min_duplication_length: args.min_length,
            max_cardinality: args.max_cardinality,

            reverse: args.reverse,
            complement: args.complement,
            skip_masked: args.skip_masked,

            compute_score: args.compute_score,
            threads_count: args.threads.unwrap_or(num_cpus::get_physical()),
            trim: args.trim.clone().map(|trim| (trim[0], trim[1])),
        },
    )?;

    let out_radix = if args.out.is_none() {
        format!(
            "{}{}{}{}{}{}.json",
            &args.prefix,
            radix,
            if args.reverse || args.complement {
                "_"
            } else {
                ""
            },
            if args.reverse { "R" } else { "" },
            if args.complement { "C" } else { "" },
            &args
                .trim
                .map(|trim| format!("_{}-{}", trim[0], trim[1]))
                .unwrap_or_default()
        )
    } else {
        args.out.unwrap()
    };

    let out_filename = asgart::utils::make_out_filename(Some(&out_radix), "", "json")
        .to_str()
        .unwrap()
        .to_owned();
    let mut out = std::fs::File::create(&out_filename)
        .with_context(|| format!("Unable to create `{}`", out_filename))?;
    let exporter = Box::new(exporters::JSONExporter) as Box<dyn exporters::Exporter>;
    exporter.save(&result, &mut out)?;
    info!(
        "{}",
        style(format!("Result written to {}", &out_filename)).bold()
    );
    Ok(())
}

fn search_duplications(strands_files: &[String], settings: RunSettings) -> Result<RunResult> {
    let total = Instant::now();

    info!("Preprocessing data");
    let (trim, to_process, mut strand) =
        prepare_data(strands_files, settings.skip_masked, settings.trim)?;

    let mut steps: Vec<Box<dyn Step>> = Vec::new();
    steps.push(Box::new(SearchDuplications::new(
        &to_process,
        trim,
        settings,
    )));
    steps.push(Box::new(FilterNs {}));
    steps.push(Box::new(ReOrder {}));
    steps.push(Box::new(ReduceOverlap {}));
    if settings.compute_score {
        steps.push(Box::new(ComputeScore {}));
    }
    steps.push(Box::new(Sort {}));

    let mut result = Vec::new();
    for (i, step) in steps.iter().enumerate() {
        info!(
            "{} {}...",
            style(format!("[{}/{}]", i + 1, steps.len())).blue().bold(),
            step.name()
        );
        result = step.run(result, &mut strand);
    }

    info!(
        "{}",
        style(format!(
            "{:?} processed in {}.",
            strands_files.join(", "),
            HumanDuration(total.elapsed())
        ))
        .green()
        .bold()
    );

    let strand = StrandResult {
        name: strand.file_names.clone(),
        length: strand.map.iter().fold(0, |ax, chr| ax + chr.length),
        map: strand.map.clone(),
    };

    Ok(RunResult {
        strand: strand.clone(),
        settings: settings,
        families: result
            .iter()
            .map(|family| {
                family
                    .iter()
                    .map(|sd| SD {
                        chr_left: strand
                            .find_chr_by_pos(sd.left)
                            .and_then(|c| Some(c.name.clone()))
                            .unwrap_or("unknown".to_string()),
                        chr_right: strand
                            .find_chr_by_pos(sd.right)
                            .and_then(|c| Some(c.name.clone()))
                            .unwrap_or("unknown".to_string()),

                        global_left_position: sd.left,
                        global_right_position: sd.right,

                        chr_left_position: sd.left
                            - strand
                                .find_chr_by_pos(sd.left)
                                .and_then(|c| Some(c.position))
                                .unwrap_or(0),
                        chr_right_position: sd.right
                            - strand
                                .find_chr_by_pos(sd.right)
                                .and_then(|c| Some(c.position))
                                .unwrap_or(0),

                        left_length: sd.left_length,
                        right_length: sd.right_length,

                        left_seq: None,
                        right_seq: None,

                        identity: sd.identity,
                        reversed: sd.reversed,
                        complemented: sd.complemented,
                    })
                    .collect::<Vec<SD>>()
            })
            .collect(),
    })
}
