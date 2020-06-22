use std::fs::File;

use anyhow::Result;
use clap::*;
use log::LevelFilter;

use asgart::exporters;
use asgart::exporters::Exporter;
use asgart::logger::Logger;
use asgart::structs::*;

fn main() -> Result<()> {
    Logger::init(LevelFilter::Info).expect("Unable to initialize logger");
    let args = App::new("ASGART slice")
        .setting(AppSettings::ColoredHelp)
        .setting(AppSettings::ColorAuto)
        .setting(AppSettings::VersionlessSubcommands)
        .setting(AppSettings::UnifiedHelpMessage)
        .version(crate_version!())
        .author(crate_authors!())
        .about("asgart-slice combines multiple ASGART JSON files into a single output file in the desired format, and features functions to filter, convert and collapse data.")
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

        .arg(Arg::with_name("max-family-members")
             .short("M")
             .long("max-family-members")
             .help("Skip families with more duplicons than specified")
             .takes_value(true))

        .arg(Arg::with_name("no-inter")
             .long("no-inter")
             .help("filters out inter-fragments duplications"))
        .arg(Arg::with_name("no-intra")
             .long("no-intra")
             .help("filters out intra-fragments duplications"))

        .arg(Arg::with_name("min-length")
             .long("min-length")
             .help("filters duplicons shorter than the given argument")
             .takes_value(true))

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
        let out_filename =
            asgart::utils::make_out_filename(args.value_of("OUTPUT"), "out", &format);
        Box::new(File::create(out_filename).unwrap())
    } else {
        Box::new(std::io::stdout())
    };

    let exporter = match &format[..] {
        "json" => Box::new(exporters::JSONExporter) as Box<dyn Exporter>,
        "gff2" => Box::new(exporters::GFF2Exporter) as Box<dyn Exporter>,
        "gff3" => Box::new(exporters::GFF3Exporter) as Box<dyn Exporter>,
        format @ _ => {
            log::warn!("Unknown output format `{}`: using json instead", format);
            Box::new(exporters::JSONExporter) as Box<dyn Exporter>
        }
    };

    let mut results = RunResult::from_files(&inputs)?;

    if args.is_present("flatten") {
        results.flatten();
    }
    if args.is_present("no-direct") {
        results.remove_direct();
    }
    if args.is_present("no-reversed") {
        results.remove_reversed();
    }
    if args.is_present("no-uncomplemented") {
        results.remove_uncomplemented();
    }
    if args.is_present("no-complemented") {
        results.remove_complemented();
    }
    if args.is_present("no-inter") {
        results.remove_inter();
    }
    if args.is_present("no-intra") {
        results.remove_intra();
    }
    if args.is_present("min-length") {
        let min_length = value_t!(args, "min-length", usize).unwrap();
        results
            .families
            .iter_mut()
            .for_each(|family| family.retain(|sd| std::cmp::min(sd.left_length, sd.right_length) >= min_length));
        results.families.retain(|f| !f.is_empty());
    }
    if args.is_present("max-family-members") {
        results
            .max_family_members(value_t!(args, "max-family-members", usize).unwrap_or(100_000_000));
    }
    if args.is_present("restrict-fragments") {
        results.restrict_fragments(&values_t!(args, "restrict-fragments", String).unwrap());
    }
    if args.is_present("exclude-fragments") {
        results.exclude_fragments(&values_t!(args, "exclude-fragments", String).unwrap());
    }

    exporter.save(&results, &mut out)?;
    Ok(())
}
