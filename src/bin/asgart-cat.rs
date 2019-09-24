#[macro_use] extern crate clap;
extern crate asgart;

use std::fs::File;
use clap::{App, Arg, AppSettings};

use asgart::separator::Separatable;
use asgart::log::LevelFilter;
use asgart::exporters::Exporter;
use asgart::exporters;
use asgart::logger::Logger;
use asgart::rayon::prelude::*;
use asgart::errors::*;
use asgart::structs::*;


fn main() {
    Logger::init(LevelFilter::Info).expect("Unable to initialize logger");
    let x = run();
    if let Err(ref e) = x {
        log::error!("{}", e);
        for e in e.iter().skip(1) { println!("{}", e); }
        if let Some(backtrace) = e.backtrace() { println!("Backtrace: {:?}", backtrace); }
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let args = App::new("ASGART concatenate")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .about("asgart-cat convert multiple input JSON files into a single one. Its main intended use is to merge together multiple files resulting from complementary runs on the same dataset, e.g. direct and palindromic duplications searches.\nIt also features some other functions to collapse, convert and filter data.")
        .arg(Arg::with_name("INPUT")
             .help("Set the input file(s) to use")
             .required(true)
             .takes_value(true)
             .min_values(1))
        .arg(Arg::with_name("format")
             .short("f")
             .long("format")
             .help("Set the desired output format")
             .possible_values(&["json", "gff2", "gff3"])
             .default_value("json"))
        .arg(Arg::with_name("OUTPUT")
             .short("o")
             .long("out")
             .help("Write the result to this file. If unspecifed, write to STDOUT.")
             .takes_value(true)
             .number_of_values(1))

        .arg(Arg::with_name("no-direct")
             .long("no-direct")
             .help("filters out direct duplications"))
        .arg(Arg::with_name("no-reversed")
             .long("no-reversed")
             .help("filters out reversed duplications"))
        .arg(Arg::with_name("no-complemented")
             .long("no-complemented")
             .help("filters out complemented duplications"))
        .arg(Arg::with_name("no-uncomplemented")
             .long("no-uncomplemented")
             .help("filters out uncomplemented duplications"))

        .arg(Arg::with_name("no-inter")
             .long("no-inter")
             .help("filters out inter-fragments duplications"))
        .arg(Arg::with_name("no-intra")
             .long("no-intra")
             .help("filters out intra-fragments duplications"))

        .arg(Arg::with_name("flatten")
             .short("F")
             .long("flatten")
             .help("Merge all the smaller-than-average-plus-one-sigma fragments into a single one (useful to deal with datasets containing large numbers of small fragments)"))

        .arg(Arg::with_name("restrict-fragments")
             .long("restrict-fragments")
             .help("ignore all fragments whose name are not in the provided list")
             .takes_value(true)
             .min_values(1))
        .arg(Arg::with_name("exclude-fragments")
             .long("exclude-fragments")
             .help("ignore all fragments whose name is in the list")
             .takes_value(true)
             .min_values(1))
        .get_matches();

    let inputs = values_t!(args, "INPUT", String).unwrap();
    let format = value_t!(args, "format", String).unwrap_or("json".to_string());
    let mut out: Box<dyn std::io::Write> = if args.is_present("OUTPUT") {
        let out_filename = asgart::utils::make_out_filename(
            args.value_of("OUTPUT"),
            "out",
            &format
        );
        Box::new(File::create(out_filename).unwrap())
    } else {
        Box::new(std::io::stdout())
    };

    let exporter = match &format[..] {
        "json"     => { Box::new(exporters::JSONExporter) as Box<dyn Exporter> }
        "gff2"     => { Box::new(exporters::GFF2Exporter) as Box<dyn Exporter> }
        "gff3"     => { Box::new(exporters::GFF3Exporter) as Box<dyn Exporter> }
        format @ _ => {
            log::warn!("Unknown output format `{}`: using json instead", format);
            Box::new(exporters::JSONExporter) as Box<dyn Exporter>
        }};

    let mut results = RunResult::from_files(&inputs)?;

    if args.is_present("no-direct") {results.remove_direct();}
    if args.is_present("no-reversed") {results.remove_reversed();}
    if args.is_present("no-uncomplemented") {results.remove_uncomplemented();}
    if args.is_present("no-complemented") {results.remove_complemented();}
    if args.is_present("no-inter") {results.remove_inter();}
    if args.is_present("no-intra") {results.remove_intra();}
    if args.is_present("restrict-fragments") {
        results.restrict_fragments(&values_t!(args, "restrict-fragments", String).unwrap());
    }
    if args.is_present("exclude-fragments") {
        results.exclude_fragments(&values_t!(args, "exclude-fragments", String).unwrap());
    }
    if args.is_present("flatten") {
        let before = results.strand.map.len();
        results.flatten();
        let after = results.strand.map.len();
        log::info!("{} fragments collapsed", (before - after).separated_string());
    }

    exporter.save(&results, &mut out)?;
    Ok(())
}
