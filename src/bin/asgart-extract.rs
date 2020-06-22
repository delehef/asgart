use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::Path;
use std::process;

use anyhow::{Context, Result};
use clap::*;
use log::*;

use asgart::logger::*;
use asgart::utils;
use asgart::exporters::{Exporter, JSONExporter};

fn read_fasta(filename: &str) -> Result<Vec<u8>> {
    let mut r = Vec::new();

    let reader = bio::io::fasta::Reader::from_file(filename)
        .with_context(|| format!("Unable to open FASTA file `{}`", filename))?;

    for record in reader.records() {
        let record = record.context(format!("Unable to parse FASTA file `{}`", filename))?;
        r.extend_from_slice(record.seq());
    }

    Ok(r)
}

fn main() -> Result<()> {
    Logger::init(LevelFilter::Info).expect("Unable to initialize logger");
    let args = App::new("ASGART extract")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .about("asgart-extract pulls out duplication families from an ASGART JSON file into a serie of FASTA files, one per family.")
        .arg(Arg::with_name("INPUT")
             .help("Set the input JSON file(s) to use")
             .required(true)
             .takes_value(true)
             .min_values(1))
        .arg(Arg::with_name("locations")
             .short("l")
             .long("locations")
             .help("Where to find the original FASTA files; might take multiple values")
             .takes_value(true))

        .arg(Arg::with_name("in-place")
             .short("I")
             .long("in-json")
             .help("write the sequences directly in the input JSON file"))
        .arg(Arg::with_name("dump")
             .short("D")
             .long("dump")
             .help("dump the sequences into a set of multiFASTA files"))
        .arg(Arg::with_name("destination")
             .short("d")
             .long("destination")
             .help("where to write the multiFASTA files")
             .takes_value(true)
             .number_of_values(1))
        .get_matches();

    if !args.is_present("in-place") && !args.is_present("dump") {
        return Err(anyhow::Error::msg(
            format!("Please specify at least one of `--in-json` or `--dump`; see --help for more details")))
    }
    let input = value_t!(args, "INPUT", String).unwrap();
    let destination = format!("{}/", value_t!(args, "destination", String).unwrap_or("./".to_string()));
    if !Path::new(&destination).is_dir() {
        return Err(anyhow::Error::msg(format!("`{}` is not a valid directory", destination)))
    }
    let locations = values_t!(args, "locations", String).unwrap_or(vec![".".to_owned()]);

    info!("Reading {}...", &input);
    let mut result = asgart::structs::RunResult::from_files(&[input.clone()])?;
    info!("Done.");

    let strands_files = result
        .strand
        .name
        .split(",")
        .map(str::trim)
        .map(|name| {
            for location in &locations {
                let path = format!("{}/{}", location, name);
                if Path::new(&path).exists() {
                    return Ok(path);
                }
            }
            Err(format!(
                "Unable to find {} in the locations provided ({})",
                name,
                locations.join(", ")
            ))
        })
        .collect::<std::result::Result<Vec<_>, _>>()
        .unwrap_or_else(|msg| {
            error!("{}", msg);
            process::exit(1)
        });

    let mut strand = Vec::new();
    for strand_file in &strands_files {
        info!("Reading {}...", strand_file);
        let new_strand = read_fasta(strand_file)
            .with_context(|| format!("Unable to read FASTA file `{}`", strand_file))?;
        strand.extend(new_strand);
        info!("Done.");
    }

    if args.is_present("in-place") {
        result.families.iter_mut().for_each(|family| {
            family.iter_mut().for_each(|mut sd| {
                let left_seq =
                    strand[sd.global_left_position..sd.global_left_position + sd.left_length].to_vec();
                let mut right_seq = strand
                    [sd.global_right_position..sd.global_right_position + sd.right_length]
                    .to_vec();
                if sd.reversed { right_seq.reverse(); }
                if sd.complemented { right_seq = utils::complemented(&right_seq); }

                sd.left_seq = Some(String::from_utf8(left_seq).unwrap());
                sd.right_seq = Some(String::from_utf8(right_seq).unwrap());
            })
        });
        JSONExporter.save(&result, &mut File::create(&input).unwrap())?
    }
    if args.is_present("dump") {
        for (i, family) in result.families.iter().enumerate() {
            let family_id = format!("family-{}", i);
            let out_file_name = format!("{}{}.fa", destination, family_id);

            for (j, sd) in family.iter().enumerate() {
                let mut file = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&out_file_name)
                    .with_context(|| format!("Unable to write results to `{}`", out_file_name))?;

                let left_seq =
                    strand[sd.global_left_position..sd.global_left_position + sd.left_length].to_vec();
                let mut right_seq = strand
                    [sd.global_right_position..sd.global_right_position + sd.right_length]
                    .to_vec();
                if sd.reversed { right_seq.reverse(); }
                if sd.complemented { right_seq = utils::complemented(&right_seq); }

                file.write_all(format!(
                    ">chr:{};start:{};end:{};family:{};duplicon:{}-1;length:{}\n",
                    sd.chr_left,
                    sd.chr_left_position,
                    sd.chr_left_position + sd.left_length,
                    i, j, sd.left_length).as_bytes())
                    .with_context(|| "Unable to write results")?;
                file.write_all(&left_seq)
                    .with_context(|| "Unable to write results")?;
                file.write_all(b"\n")
                    .with_context(|| "Unable to write results")?;
                file.write_all(
                    format!(
                        ">chr:{};start:{};end:{};family:{};duplicon:{}-2;length:{}\n",
                        sd.chr_right,
                        sd.chr_right_position,
                        sd.chr_right_position + sd.right_length,
                        i, j, sd.right_length
                    ).as_bytes(),)
                    .with_context(|| "Unable to write results")?;
                file.write_all(&right_seq)
                    .with_context(|| "Unable to write results")?;
                file.write_all(b"\n")
                    .with_context(|| "Unable to write results")?;
            }
        }
    }

    Ok(())
}
