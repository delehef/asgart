extern crate rustc_serialize;
extern crate colored;
#[macro_use]
extern crate clap;
#[macro_use]
extern crate error_chain;
extern crate asgart;


use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use rustc_serialize::json;
use clap::{App, AppSettings, ArgMatches};
use colored::Colorize;
use asgart::structs::*;
use asgart::plot::Plotter;
use asgart::plot::Settings;
use asgart::plot::chord_plot::ChordPlotter;
use asgart::plot::flat_plot::FlatPlotter;

mod errors {
    error_chain!{}
}
use errors::*;

fn read_result(file: &str) -> Result<RunResult> {
    let mut f = File::open(file).chain_err(|| format!("Unable to open {}", file))?;
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    json::decode(&s).chain_err(|| "Failed to parse JSON")
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
    };
    let plotter = FlatPlotter::new(settings);
    plotter.plot();

    Ok(())
}

fn chord(args: &ArgMatches) -> Result<()> {
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
    };
    let plotter = ChordPlotter::new(settings);
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
