extern crate rustc_serialize;
extern crate colored;
#[macro_use]
extern crate clap;
extern crate asgart;
extern crate error_chain;
extern crate regex;
extern crate rand;

use error_chain::*;
use regex::Regex;
use std::io::prelude::*;
use std::io::BufReader;
use std::fs::File;
use std::path::Path;
use rustc_serialize::json;
use clap::{App, AppSettings, ArgMatches};
use colored::Colorize;
use asgart::structs::*;
use asgart::plot::*;
use asgart::plot::chord_plot::ChordPlotter;
use asgart::plot::flat_plot::FlatPlotter;
use asgart::plot::eye_plot::EyePlotter;
use asgart::errors::*;
use std::collections::HashMap;

fn read_result(file: &str) -> Result<RunResult> {
    let mut f = File::open(file).chain_err(|| format!("Unable to open {}", file))?;
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    json::decode(&s).chain_err(|| "Failed to parse JSON")
}

fn filter_sds_genes(result: &mut RunResult, genes_families: &[Vec<Gene>], threshold: usize) {
    fn _overlap((xstart, xlen): (usize, usize), (ystart, ylen): (usize, usize)) -> bool {
        let xend = xstart + xlen;
        let yend = ystart + ylen;

        (xstart >= ystart && xstart <= yend) ||
            (ystart >= xstart && ystart <= xend)
    }

    let mut r = Vec::new();
    for sd in &result.sds{
        for family in genes_families {
            for gene in family {
                for position in &gene.positions{
                    let (start, length) = match *position {
                        GenePosition::Relative { ref chr, start, length} => {
                            let chr = result.strand1.find_chr(&chr);
                            (chr.position + start, length)
                        }
                        GenePosition::Absolute { start, length }         => { (start, length) }
                    };
                    if _overlap(sd.left_part(), (start - threshold, length + 2*threshold))
                        || _overlap(sd.right_part(), (start - threshold, length + 2*threshold))
                    {
                        r.push(sd.clone());
                    }
                }
            }
        }
    }

    result.sds = r.clone();
}


fn read_gene_file(r: &RunResult, file: &str) -> Result<Vec<Gene>> {
    let f = File::open(file).chain_err(|| format!("Unable to open {}", file))?;
    let f = BufReader::new(f);
    let mut d = HashMap::new();

    let re = Regex::new(r"(.*)\+(\d+)").unwrap();
    for (i, line) in f.
        lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.is_empty())
        .enumerate() {
            let v: Vec<&str> = line.split(';').collect();
            if v.len() != 3 { bail!("{}:L{} `{}`: incorrect format, expecting two members, found {}", file, i+1, line, v.len()); }
            let name = v[0].to_owned();

            let position = if re.is_match(v[1]) {
                let chr = re.captures(v[1]).unwrap().get(1).unwrap().as_str();
                let position = re.captures(v[1]).unwrap().get(2).unwrap().as_str().to_owned().parse::<usize>().unwrap();
                // XXX Incorrect for non-endogene runs
                if !r.strand1.has_chr(chr) {
                    bail!("Unable to find '{}'. Available: {}",
                          chr,
                          r.strand1.map.iter().fold(String::new(), |mut s, chr| { s.push_str(&format!("'{}' ", chr.name)); s })
                    );
                }
                if r.strand1.find_chr(chr).length < position {
                    bail!("{} greater than {}'s length ({})", position, chr, r.strand1.find_chr(chr).length);
                }
                GenePosition::Relative {
                    chr: chr.to_owned(),
                    start: position,
                    length: v[2].parse::<usize>().unwrap()
                }
            } else {
                GenePosition::Absolute {
                    start: v[1].parse::<usize>().unwrap(),
                    length: v[2].parse::<usize>().unwrap(),
                }
            };
            d.entry(name).or_insert_with(Vec::new).push(position);
        }

    let mut genes = Vec::new();
    for (name, positions) in &d {
        genes.push(Gene{
            name: name.to_owned(),
            positions: positions.clone(),
        })
    }

    Ok(genes)
}

fn flat(args: &ArgMatches) -> Result<()> {
    let json_file = args.value_of("FILE").unwrap();
    let result = read_result(json_file)?;


    let out_file =
        args.value_of("out")
        .and_then(|f| {
            let path = Path::new(f);
            if path.is_dir() {
                Some(path.join("out.svg").to_str().unwrap().to_owned())
            } else {
                Some(f.to_owned())
            }
        }).or(Some(format!("{}.svg", Path::new(json_file).file_stem().unwrap().to_str().unwrap()))).unwrap();

    println!("Result written to {}", out_file);
    let settings = Settings {
        result: result,
        out_file: out_file,

        size: 200.0,

        thickness: 1.0,
        color1: "#ddf983".to_owned(),
        color2: "#ddf983".to_owned(),

        min_length: if args.is_present("min_length") { value_t!(args, "min_length", usize).unwrap_or_else(|e| e.exit()) } else { 5000 },
        min_identity: 00.0, // TODO XXX
        plot_if_reversed: args.is_present("reversed"),
        plot_if_translated: args.is_present("translated"),

        gene_tracks: Vec::new(),
    };
    let plotter = FlatPlotter::new(settings);
    plotter.plot();

    Ok(())
}

fn chord(args: &ArgMatches) -> Result<()> {
    let json_file = args.value_of("FILE").unwrap();
    let mut result = read_result(json_file)?;
    let genes_tracks: Result<Vec<_>> = match args.values_of("genes") {
        Some(x) => { x
                     .map(|gene_track| read_gene_file(&result, gene_track))
                     .collect()
        }
        None    => Ok(Vec::new())
    };
    if let Err(e) = genes_tracks {return Err(e);}
    let genes_tracks = genes_tracks.unwrap();

    if args.is_present("filter_genes") {
        let threshold = value_t!(args, "filter_genes", usize).unwrap_or_else(|_| 100_000);
        filter_sds_genes(&mut result, &genes_tracks, threshold);
    }

    let out_file = args
        .value_of("out")
        .and_then(|f| {
            let path = Path::new(f);
            if path.is_dir() {
                Some(path.join("out.svg").to_str().unwrap().to_owned())
            } else {
                Some(f.to_owned())
            }
        })
        .or(Some(format!("{}.svg", Path::new(json_file).file_stem().unwrap().to_str().unwrap()))).unwrap();


    let settings = Settings {
        result: result,
        out_file: out_file,

        size: 200.0,

        thickness: 0.5,
        color1: "#ddf983".to_owned(),
        color2: "#ddf983".to_owned(),

        min_length: if args.is_present("min_length") { value_t!(args, "min_length", usize).unwrap_or_else(|e| e.exit()) } else { 5000 },
        min_identity: if args.is_present("min_identity") { value_t!(args, "min_identity", f32).unwrap_or_else(|e| e.exit()) } else { 90.0 },
        plot_if_reversed: args.is_present("reversed"),
        plot_if_translated: args.is_present("translated"),

        gene_tracks: genes_tracks,
    };
    let plotter = ChordPlotter::new(settings);
    plotter.plot();

    Ok(())
}



fn eye(args: &ArgMatches) -> Result<()> {
    let json_file = args.value_of("FILE").unwrap();
    let result = read_result(json_file)?;

    let out_file =
        args.value_of("out")
        .and_then(|f| {
            let path = Path::new(f);
            if path.is_dir() {
                Some(path.join("out.svg").to_str().unwrap().to_owned())
            } else {
                Some(f.to_owned())
            }
        }).or(Some(format!("{}.svg", Path::new(json_file).file_stem().unwrap().to_str().unwrap()))).unwrap();
    println!("Result written to {}", out_file);

    let settings = Settings {
        result: result,
        out_file: out_file,

        size: 200.0,

        thickness: 0.5,
        color1: "#ddf983".to_owned(),
        color2: "#ddf983".to_owned(),

        min_length: if args.is_present("min_length") { value_t!(args, "min_length", usize).unwrap_or_else(|e| e.exit()) } else { 5000 },
        min_identity: 90.0, // TODO
        plot_if_reversed: args.is_present("reversed"),
        plot_if_translated: args.is_present("translated"),

        gene_tracks: Vec::new(),
    };
    println!("Plotting EYE");
    let plotter = EyePlotter::new(settings);
    plotter.plot();

    Ok(())
}

fn run() -> Result<()> {
    let yaml = load_yaml!("plot.yaml");
    let args = App::from_yaml(yaml)
        .setting(AppSettings::SubcommandRequiredElseHelp)
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .get_matches();

    match args.subcommand_name() {
        Some("chord")   => chord(args.subcommand_matches("chord").unwrap()),
        Some("flat")    => flat(args.subcommand_matches("flat").unwrap()),
        Some("eye")    => eye(args.subcommand_matches("eye").unwrap()),
        None            => Ok(()),
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
            println!("backtrace: {:?}", backtrace);
        }
        std::process::exit(1);
    }
}
