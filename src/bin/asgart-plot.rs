use std::{
    collections::HashMap,
    fs::File,
    io::{prelude::*, BufReader},
    path::Path,
};

use anyhow::{anyhow, Context, Result};
use clap::*;
use regex::Regex;

use asgart::{
    plot::{
        chord_plot::ChordPlotter, circos_plot::CircosPlotter, colorizers::*,
        flat_plot::FlatPlotter, genome_plot::GenomePlotter, rosary_plot::RosaryPlotter, *,
    },
    structs::*,
};

fn filter_families_in_features(
    result: &mut RunResult,
    features_families: &[Vec<Feature>],
    threshold: usize,
) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) || (ystart >= xstart && ystart <= xend)
    }

    result.families = result
        .families
        .clone()
        .into_iter()
        .filter(|family| {
            family.iter().any(|sd| {
                features_families.iter().any(|feature_family| {
                    feature_family.iter().any(|feature| {
                        for position in &feature.positions {
                            let (start, length) = match *position {
                                FeaturePosition::Relative {
                                    ref chr,
                                    start,
                                    length,
                                } => {
                                    let chr = result
                                        .strand
                                        .find_chr(&chr)
                                        .expect(&format!("Unable to find fragment `{}`", chr));
                                    (chr.position + start, length)
                                }
                                FeaturePosition::Absolute { start, length } => (start, length),
                            };
                            if _overlap(sd.left_part(), (start - threshold, length + 2 * threshold))
                                || _overlap(
                                    sd.right_part(),
                                    (start - threshold, length + 2 * threshold),
                                )
                            {
                                return true;
                            }
                        }
                        false
                    })
                })
            })
        })
        .collect();
}

fn filter_duplicons_in_features(
    result: &mut RunResult,
    features_families: &[Vec<Feature>],
    threshold: usize,
) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) || (ystart >= xstart && ystart <= xend)
    }

    let mut families = result.families.clone();
    families.iter_mut().for_each(|family| {
        family.retain(|sd| {
            features_families.iter().any(|feature_family| {
                feature_family.iter().any(|feature| {
                    for position in &feature.positions {
                        let (start, length) = match *position {
                            FeaturePosition::Relative {
                                ref chr,
                                start,
                                length,
                            } => {
                                let chr = result
                                    .strand
                                    .find_chr(&chr)
                                    .expect(&format!("Unable to find fragment `{}`", chr));
                                (chr.position + start, length)
                            }
                            FeaturePosition::Absolute { start, length } => (start, length),
                        };
                        if _overlap(sd.left_part(), (start - threshold, length + 2 * threshold))
                            || _overlap(
                                sd.right_part(),
                                (start - threshold, length + 2 * threshold),
                            )
                        {
                            return true;
                        }
                    }
                    false
                })
            })
        })
    });
    result.families = families;
}

fn filter_features_in_sds(
    result: &mut RunResult,
    features_families: &mut [Vec<Feature>],
    threshold: usize,
) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) || (ystart >= xstart && ystart <= xend)
    }

    features_families.iter_mut().for_each(|ref mut family| {
        family.retain(|feature| {
            feature.positions.iter().any(|p| {
                let (start, length) = match *p {
                    FeaturePosition::Relative {
                        ref chr,
                        start,
                        length,
                    } => {
                        let chr = result
                            .strand
                            .find_chr(chr)
                            .unwrap_or_else(|| panic!("Unable to find fragment `{}`", chr));
                        (chr.position + start, length)
                    }
                    FeaturePosition::Absolute { start, length } => (start, length),
                };

                result.families.iter().any(|family| {
                    family.iter().any(|sd| {
                        _overlap(sd.left_part(), (start - threshold, length + 2 * threshold))
                            || _overlap(
                                sd.right_part(),
                                (start - threshold, length + 2 * threshold),
                            )
                    })
                })
            })
        })
    });
}

fn read_feature_file(r: &RunResult, file: &str) -> Result<Vec<Feature>> {
    let path = Path::new(file);
    let extension: &str = path.extension().unwrap().to_str().unwrap();

    match extension {
        "gff3" => read_gff3_feature_file(r, file),
        _ => read_custom_feature_file(r, file),
    }
}

fn read_gff3_feature_file(_r: &RunResult, file: &str) -> Result<Vec<Feature>> {
    let f = File::open(file).with_context(|| format!("Unable to open {}", file))?;
    let f = BufReader::new(f);

    let r = f
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .fold(Vec::new(), |mut ax, l| {
            let l = l.split('\t').collect::<Vec<&str>>();
            let start = l[3].parse::<usize>().unwrap();
            let end = l[4].parse::<usize>().unwrap();

            let name = if l[8].contains("Name=") {
                l[8].split(';')
                    .find(|cx| cx.contains("Name"))
                    .unwrap()
                    .split('=')
                    .nth(1)
                    .unwrap() // unwrap is safe because we check for "Name="
                    .to_string()
            } else {
                l[8].to_owned()
            };

            let feature = Feature {
                name,
                positions: vec![FeaturePosition::Relative {
                    chr: l[0].to_string(),
                    start: start,
                    length: end - start,
                }],
            };
            ax.push(feature);
            ax
        });

    Ok(r)
}

fn read_custom_feature_file(r: &RunResult, file: &str) -> Result<Vec<Feature>> {
    let f = File::open(file).with_context(|| format!("Unable to open {}", file))?;
    let f = BufReader::new(f);
    let mut d = HashMap::new();

    let re = Regex::new(r"(.*)\+(\d+)").unwrap();
    for (i, line) in f
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.is_empty())
        .filter(|l| !l.starts_with('#'))
        .enumerate()
    {
        let v: Vec<&str> = line.split(';').collect();
        if v.len() != 3 {
            return Err(anyhow!(
                "{}:L{} `{}`: incorrect format, expecting two members, found {}",
                file,
                i + 1,
                line,
                v.len()
            ));
        }
        let name = v[0].to_owned();

        let position = if re.is_match(v[1]) {
            let chr_name = re.captures(v[1]).unwrap().get(1).unwrap().as_str();
            let position = re
                .captures(v[1])
                .unwrap()
                .get(2)
                .unwrap()
                .as_str()
                .to_owned()
                .parse::<usize>()
                .unwrap();
            // XXX Incorrect for non-endofeature runs
            let chr = r
                .strand
                .find_chr(chr_name)
                .unwrap_or_else(|| panic!("Unable to find fragment `{}`", chr_name));
            if chr.length < position {
                return Err(anyhow!(
                    "{} greater than {} length ({})",
                    position,
                    chr.name,
                    chr.length
                ));
            }
            FeaturePosition::Relative {
                chr: chr.name.to_owned(),
                start: position,
                length: v[2].parse::<usize>().unwrap(),
            }
        } else {
            FeaturePosition::Absolute {
                start: v[1].parse::<usize>().unwrap(),
                length: v[2].parse::<usize>().unwrap(),
            }
        };
        d.entry(name).or_insert_with(Vec::new).push(position);
    }

    let mut features = Vec::new();
    for (name, positions) in &d {
        features.push(Feature {
            name: name.to_owned(),
            positions: positions.clone(),
        })
    }

    Ok(features)
}

#[derive(Parser)]
#[command(
    name = "ASGART slice",
    author,
    version,
    about = "Generate plots from ASGART results"
)]
struct Args {
    #[arg()]
    /// Sets the input file(s) to use. If not specified, JSON data will be expected from STDIN
    files: Option<Vec<String>>,

    #[command(flatten)]
    verbose: clap_verbosity_flag::Verbosity,

    #[arg(long)]
    /// Define a non-default output file name
    out: Option<String>,

    #[arg(long, default_value = "1000")]
    /// Filter duplicons shorter than the given value
    min_length: usize,

    #[arg(long, default_value = "0")]
    /// Filter out duplicons with a lesser identity than the given value
    min_identity: f32,

    #[arg(long, default_value = "1.0")]
    /// Filter out duplicons with a higher identity than the given value
    max_identity: f32,

    #[arg(long)]
    /// Filter out direct duplications
    no_direct: bool,

    #[arg(long)]
    /// Filter out reversed duplications
    no_reversed: bool,

    #[arg(long)]
    /// Filter out complemented duplications
    no_complemented: bool,

    #[arg(long)]
    /// Filter out non-complemented duplications
    no_uncomplemented: bool,

    #[arg(long)]
    /// Filters out inter-fragmental duplications
    no_inter: bool,

    #[arg(long)]
    /// Filters out intra-fragmental duplications
    no_intra: bool,

    #[arg(long)]
    /// Ignore all duplicons not having both arms in a fragment in the list
    restrict_fragments: Option<Vec<String>>,

    #[arg(long)]
    /// Ignore all fragments is in the given list
    exclude_fragments: Option<Vec<String>>,

    #[arg(long)]
    /// Additional feature tracks to plot
    features: Vec<String>,

    #[arg(long)]
    /// If present, do not plot duplication families further away than
    /// <filter-families> bp from features in track
    filter_families: Option<usize>,

    #[arg(long)]
    /// If present, do not plot duplicons further away than
    /// <filter-duplicons> bp from features in tracks
    filter_duplicons: Option<usize>,

    #[arg(long)]
    /// If present, do not plot features further away than
    /// <filter-features> bp from a duplicon
    filter_features: Option<usize>,

    #[arg(long, default_value = "0.1")]
    /// Plot all duplicons with at least this line thickness
    min_thickness: f64,

    #[arg(long, value_parser= ["by-type", "by-position", "by-fragment", "none"], default_value = "by-type")]
    /// Criterion on which to colorize duplicons
    colorize: String,

    #[command(subcommand)]
    /// The kind of plot to produce
    plot: Plot,
}

#[derive(Subcommand)]
enum Plot {
    /// "Plot duplications on a flat representation of the underlying fragments"
    Flat,
    /// "Plot duplications on a circo-like plot"
    Chord,
    /// "Plot duplications per chromosome on a classicaly laid out genome"
    Genome,
    /// "Generate files that can be used as input by the Circos program"
    Circos,
    /// "Plot duplications in a non-linear way for easier large-scale visualization"
    Rosary {
        #[arg(long, default_value = "0")]
        /// Two consecutive duplicons closer than the given value will be shown
        /// as a single, larger one
        clustering: usize,

        #[arg(long)]
        /// Duplications-devoid spans are represented as a string of at most
        /// 10Mbp-long beads rather than a single, larger one
        rosary: bool,
    },
}

fn main() -> Result<()> {
    let args = Args::parse();

    simple_logger::SimpleLogger::new()
        .with_level(args.verbose.log_level_filter())
        .with_colors(true)
        .init()
        .context("failed to initialize simple_logger")?;

    let (mut result, out_file) = if let Some(files) = args.files.as_ref() {
        (
            RunResult::from_files(files)?,
            asgart::utils::make_out_filename(args.out.as_deref(), &files.join("-"), ""),
        )
    } else {
        log::warn!("Reading results from STDIN");
        (
            RunResult::from_stdin()?,
            asgart::utils::make_out_filename(args.out.as_deref(), "out", ""),
        )
    };

    let mut feature_tracks = args
        .features
        .iter()
        .map(|track| read_feature_file(&result, track))
        .collect::<Result<Vec<_>>>()?;

    if args.no_direct {
        result.remove_direct();
    }
    if args.no_reversed {
        result.remove_reversed();
    }
    if args.no_uncomplemented {
        result.remove_uncomplemented();
    }
    if args.no_complemented {
        result.remove_complemented();
    }
    if args.no_inter {
        result.remove_inter();
    }
    if args.no_intra {
        result.remove_intra();
    }
    if let Some(restrict_fragments) = args.restrict_fragments.as_ref() {
        log::info!("Restricting to fragments {:?}", restrict_fragments);
        result.restrict_fragments(restrict_fragments);
    }
    if let Some(exclude_fragments) = args.exclude_fragments.as_ref() {
        log::info!("Ignoring fragments {:?}", exclude_fragments);
        result.exclude_fragments(exclude_fragments);
    }

    result.families.iter_mut().for_each(|family| {
        family.retain(|sd| sd.left_length.max(sd.right_length) >= args.min_length)
    });

    result.families.iter_mut().for_each(|family| {
        family.retain(|sd| args.min_identity <= sd.identity && sd.identity <= args.max_identity)
    });

    if let Some(filter_families) = args.filter_families {
        filter_families_in_features(&mut result, &feature_tracks, filter_families);
    }

    if let Some(filter_duplicons) = args.filter_duplicons {
        filter_duplicons_in_features(&mut result, &feature_tracks, filter_duplicons);
    }

    if let Some(filter_features) = args.filter_features {
        filter_features_in_sds(&mut result, &mut feature_tracks, filter_features);
    }

    let settings = Settings {
        out_file: out_file.to_str().unwrap().to_owned(),

        size: 200.0,
        min_thickness: args.min_thickness,
        color1: "#ff5b00".to_owned(),
        color2: "#00b2ae".to_owned(),

        feature_tracks,
    };

    let colorizer = match args.colorize.as_str() {
        "by-type" => {
            Box::new(TypeColorizer::new((1.0, 0.36, 0.0), (0.0, 0.70, 0.68))) as Box<dyn Colorizer>
        }
        "by-position" => Box::new(PositionColorizer::new(&result)) as Box<dyn Colorizer>,
        "by-fragment" => Box::new(FragmentColorizer::new(&result)) as Box<dyn Colorizer>,
        "none" => {
            Box::new(TypeColorizer::new((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))) as Box<dyn Colorizer>
        }
        _ => unreachable!(),
    };

    match args.plot {
        Plot::Flat => ChordPlotter::new(settings, result, colorizer).plot(),
        Plot::Chord => FlatPlotter::new(settings, result, colorizer).plot(),
        Plot::Genome => GenomePlotter::new(settings, result, colorizer).plot(),
        Plot::Circos => CircosPlotter::new(settings, result, colorizer).plot(),
        Plot::Rosary { clustering, rosary } => {
            RosaryPlotter::new(settings, result, colorizer, clustering, rosary).plot()
        }
    }
}
