extern crate serde_json;
extern crate colored;
#[macro_use]
extern crate clap;
extern crate asgart;
extern crate error_chain;
extern crate regex;
extern crate rand;
extern crate bio;

use error_chain::*;
use regex::Regex;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::path::Path;
use clap::{App, AppSettings};
use colored::Colorize;
use asgart::structs::*;
use asgart::plot::*;
use asgart::plot::chord_plot::ChordPlotter;
use asgart::plot::flat_plot::FlatPlotter;
use asgart::plot::circos_plot::CircosPlotter;
use asgart::plot::genome_plot::GenomePlotter;
use asgart::errors::*;
use std::collections::HashMap;

fn read_result(filename: &str) -> Result<RunResult> {
    let mut f = File::open(filename).chain_err(|| format!("Unable to open {}", filename))?;
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    serde_json::from_str(&s).chain_err(|| "Failed to parse JSON")
}

fn read_results(filenames: &[&str]) -> Result<RunResult> {
    let results = filenames
        .iter()
        .map(|filename| read_result(filename))
        .collect::<std::result::Result<Vec<RunResult>, _>>()?;

    for result in &results {
        if result.strand1.name != results[0].strand1.name || result.strand2.name != results[0].strand2.name {
            bail!(format!("trying to combine ASGART files containing different karyotypes: `{}/{}` and `{}/{}`",
                          result.strand1.name, result.strand2.name,
                          results[0].strand1.name, results[0].strand2.name
            ));
        }
    }

    Ok(RunResult {
        settings: results[0].settings.clone(),
        strand1:  results[0].strand1.clone(),
        strand2:  results[0].strand2.clone(),
        families: results.iter().flat_map(|ref r| r.families.iter()).cloned().collect::<Vec<_>>(),
    })
}

fn filter_sds_in_features(result: &mut RunResult, features_families: &[Vec<Feature>], threshold: usize) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) ||
            (ystart >= xstart && ystart <= xend)
    }

    result.families = result.families.clone().into_iter()
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
                                                let chr = result.strand1.find_chr(&chr);
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
                                      let chr = result.strand1.find_chr(&chr);
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

            let name = if l[8].contains("Name") {
                l[8]
                    .split(';')
                    .find(|cx| cx.contains("Name")).unwrap()
                    .split('=')
                    .nth(1).unwrap()
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
                let chr = re.captures(v[1]).unwrap().get(1).unwrap().as_str();
                let position = re.captures(v[1]).unwrap().get(2).unwrap().as_str().to_owned().parse::<usize>().unwrap();
                // XXX Incorrect for non-endofeature runs
                if !r.strand1.has_chr(chr) {
                    bail!("Unable to find '{}'. Available: {}",
                          chr,
                          r.strand1.map.iter().fold(String::new(), |mut s, chr| { s.push_str(&format!("'{}' ", chr.name)); s })
                    );
                }
                if r.strand1.find_chr(chr).length < position {
                    bail!("{} greater than {}'s length ({})", position, chr, r.strand1.find_chr(chr).length);
                }
                FeaturePosition::Relative {
                    chr: chr.to_owned(),
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
    let yaml = load_yaml!("plot.yaml");
    let args = App::from_yaml(yaml)
        .version(crate_version!())
        .author(crate_authors!())
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .get_matches();

    let json_files = args.values_of("FILE").unwrap().collect::<Vec<_>>();
    let mut result = read_results(&json_files)?;

    let out_file = args
        .value_of("out")
        .and_then(|f| {
            let path = Path::new(f);
            if path.is_dir() {
                Some(path.join("out").to_str().unwrap().to_owned())
            } else {
                Some(f.to_owned())
            }
        })
        .or(Some(Path::new(json_files[0]).file_stem().unwrap().to_str().unwrap().to_owned()))
        .unwrap();

    let features_tracks: Result<Vec<_>> = match args.values_of("features") {
        Some(x) => { x
                     .map(|feature_track| read_feature_file(&result, feature_track))
                     .collect()
        }
        None    => Ok(Vec::new())
    };
    if let Err(e) = features_tracks {return Err(e);}
    let mut features_tracks = features_tracks.unwrap();


    if args.is_present("no-direct") {
        result.families.iter_mut().for_each(|family| family.retain(|sd| sd.reversed))
    }
    if args.is_present("no-reversed") {
        result.families.iter_mut().for_each(|family| family.retain(|sd| !sd.reversed))
    }
    if args.is_present("no-uncomplemented") {
        result.families.iter_mut().for_each(|family| family.retain(|sd| sd.complemented))
    }
    if args.is_present("no-complemented") {
        result.families.iter_mut().for_each(|family| family.retain(|sd| !sd.complemented))
    }

    let min_length   = value_t!(args, "min_length", usize).unwrap();
    result.families.iter_mut().for_each(|family| family.retain(|sd| sd.length >= min_length));
    let min_identity = value_t!(args, "min_identity", f32).unwrap();
    result.families.iter_mut().for_each(|family| family.retain(|sd| sd.identity >= min_identity));

    if args.is_present("filter_duplications") {
        filter_sds_in_features(&mut result, &features_tracks, value_t!(args, "filter_duplications", usize).unwrap());
    }
    if args.is_present("filter_features") {
        filter_features_in_sds(&mut result, &mut features_tracks, value_t!(args, "filter_features", usize).unwrap());
    }

    let settings = Settings {
        out_file:                out_file,

        size:                    200.0,
        min_thickness:           0.1,
        color1:                  "#ff5b00".to_owned(),
        color2:                  "#00b2ae".to_owned(),

        feature_tracks:          features_tracks,
    };


    match args.value_of("PLOT-TYPE") {
        Some("chord")   => ChordPlotter::new(settings, result).plot(),
        Some("flat")    => FlatPlotter::new(settings, result).plot(),
        Some("genome")  => GenomePlotter::new(settings, result).plot(),
        Some("circos")  => CircosPlotter::new(settings, result).plot(),
        // Some("eye")    => eye(args.subcommand_matches("eye").unwrap()),
        _               => unreachable!(),
    }
}

fn main() {
    if let Err(ref e) = run() {
        println!("{} {}", "Error: ".red(), e);
        for e in e.iter().skip(1) {
            println!("{}", e);
        }
        if let Some(backtrace) = e.backtrace() {
            println!("Backtrace: {:?}", backtrace);
        }
        std::process::exit(1);
    }
}
