#[macro_use]
extern crate log;
#[macro_use]
extern crate clap;
extern crate asgart;

use std::fs::OpenOptions;
use std::io::Write;
use std::path::Path;
use std::process;

use clap::{App, AppSettings, Arg};

use asgart::log::LevelFilter;
use asgart::logger::Logger;
use asgart::errors::*;
use asgart::utils;


fn read_fasta(filename: &str) -> Result<(Vec<u8>)> {
    let mut r = Vec::new();

    let reader = bio::io::fasta::Reader::from_file(filename)
        .chain_err(|| format!("Unable to open `{}`", filename))?;

    for record in reader.records() {
        let record = record.chain_err(|| {
            format!("Unable to read {:?}: not a FASTA file", Path::new(filename).file_name().unwrap())
        })?;

        r.extend_from_slice(record.seq());
    }

    Ok(r)
}

fn main() {
    Logger::init(LevelFilter::Info).expect("Unable to initialize logger");
    let x = run();
    if let Err(ref e) = x {
        error!("{}", e);
        for e in e.iter().skip(1) {
            println!("{}", e);
        }
        if let Some(backtrace) = e.backtrace() {
            println!("Backtrace: {:?}", backtrace);
        }
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let args = App::new("ASGART extract")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .about("ASGART extract pulls out duplication families as described in an ASGART JSON file into a serie of FASTA files, one per family.")
        .arg(Arg::with_name("INPUT")
             .help("Set the input JSON file(s) to use")
             .required(true)
             .takes_value(true)
             .min_values(1))
        .arg(Arg::with_name("prefix")
             .short("p")
             .long("prefix")
             .help("A prefix to append to the output files")
             .takes_value(true)
             .number_of_values(1))
        .arg(Arg::with_name("locations")
             .short("l")
             .long("locations")
             .help("Where to find the original FASTA files; might take multiple values")
             .takes_value(true))
        .get_matches();

    let input = value_t!(args, "INPUT", String).unwrap();
    let prefix = value_t!(args, "prefix", String).unwrap_or("".to_string());
    let locations = values_t!(args, "locations", String).unwrap_or(vec!["".to_owned()]);

    info!("Reading {}...", &input);
    let result = asgart::structs::RunResult::from_files(&[input])?;
    info!("Done.");

    let strands_files = result.strand.name.split(",")
        .map(str::trim)
        .map(|name| {
            for location in &locations {
                let path = format!("{}/{}", location, name);
                if Path::new(&path).exists() {
                    return Ok(path)
                }
            }
            Err(format!("Unable to find {} in the locations provided ({})", name, locations.join(", ")))
        })
        .collect::<std::result::Result<Vec<_>, _>>()
        .unwrap_or_else(|msg| {error!("{}", msg); process::exit(1); Vec::new()});

    let mut strand = Vec::new();
    for strand_file in &strands_files {
        info!("Reading {}...", strand_file);
        let new_strand = read_fasta(strand_file).chain_err(|| format!("Unable to parse `{}`", strand_file))?;
        strand.extend(new_strand);
        info!("Done.");
    }

    info!("Writing results...");
    for (i, family) in result.families.iter().enumerate() {
        let family_id = format!("family-{}", i);
        let out_file_name = format!("{}{}.fa", prefix, family_id);

        for (j, sd) in family.iter().enumerate() {
            let mut file = OpenOptions::new().append(true).create(true).open(&out_file_name)
                .chain_err(|| format!("Unable to write in `{}`", out_file_name))?;

            let left_seq = strand[sd.global_left_position .. sd.global_left_position + sd.left_length].to_vec();
            let mut right_seq = strand[sd.global_right_position .. sd.global_right_position + sd.right_length].to_vec();
            if sd.reversed {right_seq.reverse();}
            if sd.complemented {right_seq = utils::complemented(&right_seq);}

            file.write_all(format!(">{} {}-{} duplicon-{}-1\n",
                                   sd.chr_left, sd.chr_left_position, sd.chr_left_position+sd.left_length , j).as_bytes())
                .chain_err(|| "Unable to write in file.")?;
            file.write_all(&left_seq).chain_err(|| "Unable to write in file.")?;
            file.write_all(b"\n").chain_err(|| "Unable to write in file.")?;
            file.write_all(format!(">{} {}-{} duplicon-{}-2\n",
                                   sd.chr_right, sd.chr_right_position, sd.chr_right_position+sd.right_length, j).as_bytes())
                .chain_err(|| "Unable to write in file.")?;
            file.write_all(&right_seq).chain_err(|| "Unable to write in file.")?;
            file.write_all(b"\n").chain_err(|| "Unable to write in file.")?;
        }
    }
    info!("Done.");


    Ok(())
}
