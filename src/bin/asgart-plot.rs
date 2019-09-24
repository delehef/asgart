#[macro_use] extern crate clap;
extern crate serde_json;
extern crate colored;
extern crate asgart;
extern crate error_chain;
extern crate regex;
extern crate rand;
extern crate bio;

use error_chain::*;
use regex::Regex;
use std::io::BufReader;
use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use std::collections::HashMap;

use clap::{App, AppSettings};
use colored::Colorize;

use asgart::structs::*;
use asgart::plot::*;
use asgart::plot::chord_plot::ChordPlotter;
use asgart::plot::flat_plot::FlatPlotter;
use asgart::plot::circos_plot::CircosPlotter;
use asgart::plot::genome_plot::GenomePlotter;
use asgart::plot::colorizers::*;
use asgart::errors::*;
use asgart::logger::Logger;
use asgart::log::LevelFilter;


fn filter_families_in_features(result: &mut RunResult, features_families: &[Vec<Feature>], threshold: usize) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) ||
            (ystart >= xstart && ystart <= xend)
    }

    result.families = result.families
        .clone()
        .into_iter()
        .filter(
            |family|
            family.iter()
                .any(|sd| {
                    features_families.iter()
                        .any(|feature_family| {
                            feature_family.iter()
                                .any(|feature| {
                                    for position in &feature.positions{
                                        let (start, length) = match *position {
                                            FeaturePosition::Relative { ref chr, start, length} => {
                                                let chr = result.strand.find_chr(&chr).expect(&format!("Unable to find fragment `{}`", chr));
                                                (chr.position + start, length)
                                            }
                                            FeaturePosition::Absolute { start, length }         => {
                                                (start, length)
                                            }
                                        };
                                        if _overlap(sd.left_part(), (start - threshold, length + 2*threshold))
                                            || _overlap(sd.right_part(), (start - threshold, length + 2*threshold))
                                        {
                                            return true;
                                        }
                                    }
                                    false
                                })
                        })
                })
        )
        .collect();
}


fn filter_duplicons_in_features(result: &mut RunResult, features_families: &[Vec<Feature>], threshold: usize) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) ||
            (ystart >= xstart && ystart <= xend)
    }

    let mut families = result.families.clone();
    families
        .iter_mut()
        .for_each(
            |family|
            family
                .retain(|sd| {
                    features_families.iter()
                        .any(|feature_family| {
                            feature_family.iter()
                                .any(|feature| {
                                    for position in &feature.positions{
                                        let (start, length) = match *position {
                                            FeaturePosition::Relative { ref chr, start, length} => {
                                                let chr = result.strand.find_chr(&chr).expect(&format!("Unable to find fragment `{}`", chr));
                                                (chr.position + start, length)
                                            }
                                            FeaturePosition::Absolute { start, length }         => {
                                                (start, length)
                                            }
                                        };
                                        if _overlap(sd.left_part(), (start - threshold, length + 2*threshold))
                                            || _overlap(sd.right_part(), (start - threshold, length + 2*threshold))
                                        {
                                            return true;
                                        }
                                    }
                                    false
                                })
                        })
                })
        );
    result.families = families;
}

fn filter_features_in_sds(result: &mut RunResult, features_families: &mut Vec<Vec<Feature>>, threshold: usize) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) ||
            (ystart >= xstart && ystart <= xend)
    }

    features_families
        .iter_mut()
        .for_each(|ref mut family|
                  family
                  .retain(|feature| {
                      feature.positions
                          .iter()
                          .any(|p| {
                              let (start, length) = match *p {
                                  FeaturePosition::Relative { ref chr, start, length} => {
                                      let chr = result.strand.find_chr(&chr).expect(&format!("Unable to find fragment `{}`", chr));
                                      (chr.position + start, length)
                                  }
                                  FeaturePosition::Absolute { start, length }         => { (start, length) }
                              };

                              result.families.iter().any(
                                  |family| family.iter().any(
                                      |ref sd|
                                      _overlap(sd.left_part(), (start - threshold, length + 2*threshold))
                                          || _overlap(sd.right_part(), (start - threshold, length + 2*threshold))
                                  )
                              )
                          })
                  })
        );
}

fn read_feature_file(r: &RunResult, file: &str) -> Result<Vec<Feature>> {
    let path = Path::new(file);
    let extension: &str = path.extension().unwrap().to_str().unwrap();

    match extension {
        "gff3" => {read_gff3_feature_file(r, file)}
        _      => {read_custom_feature_file(r, file)}
    }
}

fn read_gff3_feature_file(_r: &RunResult, file: &str) -> Result<Vec<Feature>> {
    let f = File::open(file).chain_err(|| format!("Unable to open {}", file))?;
    let f = BufReader::new(f);

    let r = f.lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.is_empty() && !l.starts_with('#'))
        .fold(Vec::new(), |mut ax, l| {
            let l = l.split('\t').collect::<Vec<&str>>();
            let start = l[3].parse::<usize>().unwrap();
            let end = l[4].parse::<usize>().unwrap();

            let name = if l[8].contains("Name=") {
                l[8]
                    .split(';')
                    .find(|cx| cx.contains("Name")).unwrap()
                    .split('=')
                    .nth(1).unwrap() // unwrap is safe because we check for "Name="
                    .to_string()
            } else {
                l[8].to_owned()
            };

            let feature = Feature {
                name: name,
                positions: vec![FeaturePosition::Relative{
                    chr: l[0].to_string(),
                    start:  start,
                    length: end - start,
                }],
            } ;
            ax.push(feature);
            ax
        });

    Ok(r)
}

fn read_custom_feature_file(r: &RunResult, file: &str) -> Result<Vec<Feature>> {
    let f = File::open(file).chain_err(|| format!("Unable to open {}", file))?;
    let f = BufReader::new(f);
    let mut d = HashMap::new();

    let re = Regex::new(r"(.*)\+(\d+)").unwrap();
    for (i, line) in f.
        lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.is_empty())
        .filter(|l| !l.starts_with('#'))
        .enumerate() {
            let v: Vec<&str> = line.split(';').collect();
            if v.len() != 3 { bail!("{}:L{} `{}`: incorrect format, expecting two members, found {}", file, i+1, line, v.len()); }
            let name = v[0].to_owned();

            let position = if re.is_match(v[1]) {
                let chr_name = re.captures(v[1]).unwrap().get(1).unwrap().as_str();
                let position = re.captures(v[1]).unwrap().get(2).unwrap().as_str().to_owned().parse::<usize>().unwrap();
                // XXX Incorrect for non-endofeature runs
                let chr = r.strand.find_chr(chr_name).expect(&format!("Unable to find fragment `{}`", chr_name));
                if chr.length < position {
                    bail!("{} greater than {} length ({})", position, chr.name, chr.length);
                }
                FeaturePosition::Relative {
                    chr: chr.name.to_owned(),
                    start: position,
                    length: v[2].parse::<usize>().unwrap()
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
        features.push(Feature{
            name: name.to_owned(),
            positions: positions.clone(),
        })
    }

    Ok(features)
}

fn run() -> Result<()> {
    let yaml = load_yaml!("asgart-plot.yaml");
    let args = App::from_yaml(yaml)
        .version(crate_version!())
        .author(crate_authors!())
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .get_matches();
    Logger::init(match args.occurrences_of("verbose") {
        0 => { LevelFilter::Info }
        1 => { LevelFilter::Debug }
        2 => { LevelFilter::Trace }
        _ => { LevelFilter::Trace }
    }).chain_err(|| "Unable to initialize logger")?;

    let json_files = values_t!(args, "FILE", String).unwrap();
    let mut result = RunResult::from_files(&json_files)?;

    let out_file = asgart::utils::make_out_filename(args.value_of("out"), &json_files.join("-"), "");

    let features_tracks: Result<Vec<_>> = match args.values_of("features") {
        Some(x) => { x
                    .map(|feature_track| read_feature_file(&result, feature_track))
                    .collect()
        }
        None    => Ok(Vec::new())
    };
    if let Err(e) = features_tracks {return Err(e);}
    let mut features_tracks = features_tracks.unwrap();


    if args.is_present("no-direct") {result.remove_direct();}
    if args.is_present("no-reversed") {result.remove_reversed();}
    if args.is_present("no-uncomplemented") {result.remove_uncomplemented();}
    if args.is_present("no-complemented") {result.remove_complemented();}
    if args.is_present("no-inter") {result.remove_inter();}
    if args.is_present("no-intra") {result.remove_intra();}
    if args.is_present("restrict-fragments") {
        let to_keep = values_t!(args, "restrict-fragments", String).unwrap();
        log::info!("Restricting to fragments {:?}", &to_keep);
        result.restrict_fragments(&to_keep);
    }
    if args.is_present("exclude-fragments") {
        let to_remove = &values_t!(args, "exclude-fragments", String).unwrap();
        log::info!("Ignoring fragments {:?}", &to_remove);
        result.exclude_fragments(&to_remove);
    }


    let min_length   = value_t!(args, "min_length", usize).unwrap();
    result.families.iter_mut().for_each(|family| family.retain(|sd| std::cmp::max(sd.left_length, sd.right_length) >= min_length));
    let min_identity = value_t!(args, "min_identity", f32).unwrap();
    result.families.iter_mut().for_each(|family| family.retain(|sd| sd.identity >= min_identity));

    if args.is_present("filter_families") {
        filter_families_in_features(
            &mut result, &features_tracks,
            value_t!(args, "filter_families", usize).unwrap()
        );
    }
    if args.is_present("filter_duplicons") {
        filter_duplicons_in_features(
            &mut result, &features_tracks,
            value_t!(args, "filter_duplicons", usize).unwrap()
        );
    }

    if args.is_present("filter_features") {
        filter_features_in_sds(
            &mut result, &mut features_tracks,
            value_t!(args, "filter_features", usize).unwrap()
        );
    }

    let settings = Settings {
        out_file:                out_file.to_str().unwrap().to_owned(),

        size:                    200.0,
        min_thickness:           value_t!(args, "min_thickness", f64).unwrap(),
        color1:                  "#ff5b00".to_owned(),
        color2:                  "#00b2ae".to_owned(),

        feature_tracks:          features_tracks,
    };


    let colorizer = match args.value_of("colorize") {
        Some("by-type")     => Box::new(TypeColorizer::new((1.0, 0.36, 0.0), (0.0, 0.70, 0.68))) as Box<dyn Colorizer>,
        Some("by-position") => Box::new(PositionColorizer::new(&result)) as Box<dyn Colorizer>,
        Some("by-fragment") => Box::new(FragmentColorizer::new(&result)) as Box<dyn Colorizer>,
        Some("none")        => Box::new(TypeColorizer::new((0.5, 0.5, 0.5), (0.5, 0.5, 0.5))) as Box<dyn Colorizer>,
        _                   => unreachable!(),
    } ;

    let r = match args.value_of("PLOT-TYPE") {
        Some("chord")   => ChordPlotter::new(settings, result, colorizer).plot(),
        Some("flat")    => FlatPlotter::new(settings, result, colorizer).plot(),
        Some("genome")  => GenomePlotter::new(settings, result, colorizer).plot(),
        Some("circos")  => CircosPlotter::new(settings, result, colorizer).plot(),
        // Some("eye")    => eye(args.subcommand_matches("eye").unwrap()),
        _               => unreachable!(),
    };
    r
}

fn main() {
    if let Err(ref e) = run() {
        println!("{} {}", "Error: ".red(), e);
        for e in e.iter().skip(1) { println!("{}", e); }
        if let Some(backtrace) = e.backtrace() {
            println!("Backtrace: {:?}", backtrace);
        }
        std::process::exit(1);
    }
}
