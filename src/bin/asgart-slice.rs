use std::fs::File;

use anyhow::{Context, Result};
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
             .min_values(1)
             .takes_value(true))
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
             .takes_value(true)
             .number_of_values(1))

        .arg(Arg::with_name("no-inter")
             .long("no-inter")
             .help("filters out inter-fragments duplications")
             .conflicts_with("no-inter-relaxed"))
        .arg(Arg::with_name("no-inter-relaxed")
             .long("no-inter-relaxed")
             .help("filters out inter-fragments duplications, except in they stand into the collapsed pseudo-chromosome")
             .conflicts_with("no-inter"))
        .arg(Arg::with_name("no-intra")
             .long("no-intra")
             .help("filters out intra-fragments duplications"))

        .arg(Arg::with_name("min-length")
             .long("min-length")
             .help("filters duplicons shorter than the given argument")
             .takes_value(true)
             .number_of_values(1))

        .arg(Arg::with_name("collapse")
             .short("C")
             .long("collapse")
             .help("Merge all the smaller-than-average-plus-one-sigma fragments into a single one (useful to deal with datasets containing large numbers of small fragments)"))

        .arg(Arg::with_name("keep-fragments")
             .long("keep-fragments")
             .help("ignore all fragments not in the list [comma-separated]")
             .takes_value(true)
             .min_values(1)
             .require_delimiter(true))
        .arg(Arg::with_name("exclude-fragments")
             .long("exclude-fragments")
             .help("ignore all fragments is in the list [comma-separated]")
             .takes_value(true)
             .min_values(1)
             .require_delimiter(true))
        .arg(Arg::with_name("regexp")
             .short("E")
             .long("--regexp")
             .help("Use regexp matching instead of literal for keep- and exclude-fragments"))
        .get_matches();

    let mut results = if args.is_present("INPUT") {
        let inputs = values_t!(args, "INPUT", String).unwrap();
        RunResult::from_files(&inputs)?
    } else {
        log::warn!("Reading results from STDIN");
        RunResult::from_stdin()?
    };

    let format = value_t!(args, "format", String).unwrap_or("json".to_string());
    let mut out: Box<dyn std::io::Write> = if args.is_present("OUTPUT") {
        let out_filename =
            asgart::utils::make_out_filename(args.value_of("OUTPUT"), "out", &format);
        Box::new(File::create(out_filename)?)
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

    if args.is_present("collapse") {
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
    if args.is_present("no-inter-relaxed") {
        results.remove_inter_relaxed();
    }
    if args.is_present("no-intra") {
        results.remove_intra();
    }
    if args.is_present("min-length") {
        let min_length = value_t!(args, "min-length", usize)?;
        results.families.iter_mut().for_each(|family| {
            family.retain(|sd| std::cmp::min(sd.left_length, sd.right_length) >= min_length)
        });
        results.families.retain(|f| !f.is_empty());
    }
    if args.is_present("max-family-members") {
        results
            .max_family_members(value_t!(args, "max-family-members", usize).unwrap_or(100_000_000));
    }
    if args.is_present("keep-fragments") {
        if args.is_present("regexp") {
            let regexp = value_t!(args, "keep-fragments", String)?;
            results
                .restrict_fragments_regexp::<&str>(&regexp)
                .with_context(|| format!("Error while compiling `{}`", regexp))?;
        } else {
            results.restrict_fragments(&values_t!(args, "keep-fragments", String).unwrap());
        }
    }
    if args.is_present("exclude-fragments") {
        if args.is_present("regexp") {
            let regexp = value_t!(args, "exclude-fragments", String)?;
            results
                .exclude_fragments_regexp::<&str>(&regexp)
                .with_context(|| format!("Error while compiling `{}`", regexp))?;
        } else {
            results.exclude_fragments(&values_t!(args, "exclude-fragments", String)?);
        }
    }

    exporter.save(&results, &mut out)?;
    Ok(())
}
