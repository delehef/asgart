extern crate rustc_serialize;
extern crate colored;
#[macro_use]
extern crate clap;
extern crate asgart;


use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use rustc_serialize::json;
use clap::{App, AppSettings, ArgMatches};
use asgart::structs::*;
use asgart::plot::Plotter;
use asgart::plot::Settings;
use asgart::plot::chord_plot::ChordPlotter;
use asgart::plot::flat_plot::FlatPlotter;


fn flat(args: &ArgMatches) {
    println!("Flat plotting {}", args.value_of("FILE").unwrap());

    let mut f = File::open(args.value_of("FILE").unwrap()).unwrap();
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    let result: RunResult = json::decode(&s).expect("Unable to parse JSON file");

    let settings = Settings {
        result: result,
        out_file: args.value_of("out").unwrap_or_else(|| "out.svg").to_owned(),

        size: 200.0,

        thickness: 1.0,
        color1: "#ddf983".to_owned(),
        color2: "#ddf983".to_owned(),

        min_length: 5000,
    };
    let plotter = FlatPlotter::new(settings);
    plotter.plot();
}

fn chord(args: &ArgMatches) {
    let json_file = args.value_of("FILE").unwrap();
    println!("Chord plotting {}", json_file);

    let mut f = File::open(json_file).expect(&format!("Unable to open {}", json_file));
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    let result: RunResult = json::decode(&s).expect("Unable to parse JSON file");


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

    println!("Writing to {}", out_file);

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
}

fn main() {
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
        Some("dotplot") => println!("'git add' was used"),
        None            => println!("No subcommand was used"),
        _               => unreachable!(),
    }
}
