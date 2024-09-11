use std::{
    fs::{File, OpenOptions},
    io::Write,
    path::Path,
    process,
};

use anyhow::{Context, Result};
use clap::*;
use log::*;

use asgart::{
    exporters::{Exporter, JSONExporter},
    utils,
};

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

#[derive(Parser)]
#[command(
    name = "ASGART extract",
    version,
    author,
    about = "asgart-extract pulls out duplication families from an ASGART JSON file into a serie of FASTA files, one per family."
)]
struct Args {
    #[arg()]
    /// The JSON file to process
    input: String,

    #[arg(short = 'l', long)]
    /// Where to find the original FASTA files; multiple values can be given
    locations: Option<Vec<String>>,

    #[arg(short = 'I', long)]
    /// Write the sequences directly into the input JSON files
    in_place: bool,

    #[arg(short = 'D', long)]
    /// Dump the sequences into multiFASTA files
    dump: bool,

    #[arg(short = 'd', long)]
    /// Where to write the output multiFASTA files
    destination: Option<String>,
}

fn main() -> Result<()> {
    simple_logger::SimpleLogger::new()
        .with_level(LevelFilter::Info)
        .with_colors(true)
        .init()
        .context("failed to initialize simple_logger")?;

    let args = Args::parse();

    if !args.in_place && !args.dump {
        return Err(anyhow::Error::msg(format!(
            "Please specify at least one of `--in-place` or `--dump`; see --help for more details"
        )));
    }
    let destination = format!("{}/", args.destination.unwrap_or("./".to_string()));
    if !Path::new(&destination).is_dir() {
        return Err(anyhow::Error::msg(format!(
            "`{}` is not a valid directory",
            destination
        )));
    }
    let locations = args.locations.unwrap_or(vec![".".to_owned()]);

    info!("Reading {}...", &args.input);
    let mut result = asgart::structs::RunResult::from_files(&[args.input.clone()])?;
    info!("Done.");

    let strands_files = result
        .strand
        .name
        .split(',')
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

    if args.in_place {
        result.families.iter_mut().for_each(|family| {
            family.iter_mut().for_each(|sd| {
                let left_seq = strand
                    [sd.global_left_position..sd.global_left_position + sd.left_length]
                    .to_vec();
                let mut right_seq = strand
                    [sd.global_right_position..sd.global_right_position + sd.right_length]
                    .to_vec();
                if sd.reversed {
                    right_seq.reverse();
                }
                if sd.complemented {
                    right_seq = utils::complemented(&right_seq);
                }

                sd.left_seq = Some(String::from_utf8(left_seq).unwrap());
                sd.right_seq = Some(String::from_utf8(right_seq).unwrap());
            })
        });
        JSONExporter.save(&result, &mut File::create(&args.input).unwrap())?
    }
    if args.dump {
        for (i, family) in result.families.iter().enumerate() {
            let family_id = format!("family-{}", i);
            let out_file_name = format!("{}{}.fa", destination, family_id);

            for (j, sd) in family.iter().enumerate() {
                let mut file = OpenOptions::new()
                    .append(true)
                    .create(true)
                    .open(&out_file_name)
                    .with_context(|| format!("Unable to write results to `{}`", out_file_name))?;

                let left_seq = strand
                    [sd.global_left_position..sd.global_left_position + sd.left_length]
                    .to_vec();
                let mut right_seq = strand
                    [sd.global_right_position..sd.global_right_position + sd.right_length]
                    .to_vec();
                if sd.reversed {
                    right_seq.reverse();
                }
                if sd.complemented {
                    right_seq = utils::complemented(&right_seq);
                }

                file.write_all(
                    format!(
                        ">chr:{};start:{};end:{};family:{};duplicon:{}-1;length:{}\n",
                        sd.chr_left,
                        sd.chr_left_position,
                        sd.chr_left_position + sd.left_length,
                        i,
                        j,
                        sd.left_length
                    )
                    .as_bytes(),
                )
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
                        i,
                        j,
                        sd.right_length
                    )
                    .as_bytes(),
                )
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
