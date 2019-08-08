#[macro_use] extern crate log;
#[macro_use] extern crate clap;
#[macro_use] extern crate asgart;

use std::fs::File;

use clap::{App, Arg, AppSettings};
use asgart::exporters::Exporter;
use asgart::exporters;
use asgart::logger::Logger;
use asgart::log::LevelFilter;

use asgart::error_chain::*;
use asgart::errors::*;

fn main() {
    Logger::init(LevelFilter::Info).expect("Unable to initialize logger");
    let x = run();
    if let Err(ref e) = x {
        error!("{}", e);
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
        .about("ASGART concatenate convert multiple input JSON files into a single one. Its main intended use is to merge together multiple files resulting from complementary runs on the same dataset, e.g. direct and palindromic duplications searches.")
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
        .get_matches();

    let inputs = values_t!(args, "INPUT", String).unwrap();
    let format = value_t!(args, "format", String).unwrap_or("json".to_string());
    let mut out: Box<dyn std::io::Write> = if args.is_present("OUTPUT") {
        let out_filename = asgart::utils::make_out_filename(args.value_of("OUTPUT"), "out", &format);
        Box::new(File::create(out_filename).expect(""))
    } else {
        Box::new(std::io::stdout())
    };
    let exporter = match &format[..] {
        "json"     => { Box::new(exporters::JSONExporter) as Box<dyn Exporter> }
        "gff2"     => { Box::new(exporters::GFF2Exporter) as Box<dyn Exporter> }
        "gff3"     => { Box::new(exporters::GFF3Exporter) as Box<dyn Exporter> }
        format @ _ => {
            warn!("Unknown output format `{}`: using json instead", format);
            Box::new(exporters::JSONExporter) as Box<dyn Exporter>
        }};

    let results = asgart::structs::RunResult::from_files(&inputs)?;
    exporter.save(&results, &mut out)?;
    Ok(())
}
