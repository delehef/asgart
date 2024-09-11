use std::fs::File;

use anyhow::{Context, Result};
use clap::*;
use log::LevelFilter;

use asgart::{
    exporters::{self, Exporter},
    structs::*,
};

#[derive(Parser)]
#[command(
    name = "ASGART slice",
    version,
    author,
    about = "asgart-slice combines multiple ASGART JSON files into a single output file in the desired format, and features functions to filter, convert and collapse data."
)]
struct Args {
    #[arg()]
    /// The input file(s) to slice
    inputs: Vec<String>,

    #[arg(short='f', long, value_parser=["json", "gff2", "gff3"], default_value="json")]
    /// Set the desired output format
    format: String,

    #[arg(short = 'o', long)]
    /// If specified, write the result to this file; otherwise, write to STDOUT
    output: Option<String>,

    #[arg(long)]
    /// Filter out direct duplications
    no_direct: bool,

    #[arg(long)]
    /// Filter out reversed duplications
    no_reversed: bool,

    #[arg(long)]
    /// Filter out complemented duplications
    no_complemented: bool,

    #[arg(long)]
    /// Filter out non-complemented duplications
    no_uncomplemented: bool,

    #[arg(short = 'M', long)]
    /// Skip families with more duplicons than specified
    max_family_members: Option<usize>,

    #[arg(long)]
    /// Filters out inter-fragmental duplications
    no_inter: bool,

    #[arg(long, conflicts_with("no_inter"))]
    /// Filters out inter-fragmental duplications, except when they lay in the
    /// collapsed pseudo-chromosome
    no_inter_relaxed: bool,

    #[arg(long)]
    /// Filters out intra-fragmental duplications
    no_intra: bool,

    #[arg(long)]
    /// Filter duplicons shorter than the given value
    min_length: Option<usize>,

    #[arg(short = 'C', long)]
    /// Merge all the smaller-than-average-plus-one-sigma fragments into a
    /// single one (useful to deal with datasets containing large numbers of
    /// small fragments)
    collapse: bool,

    #[arg(long)]
    /// Ignore all duplicons not having at least an arm in a fragment in the
    /// given list
    keep_fragments: Option<Vec<String>>,

    #[arg(long)]
    /// Ignore all duplicons not having both arms in a fragment in the list
    restrict_fragments: Option<Vec<String>>,

    #[arg(long)]
    /// Ignore all fragments is in the given list
    exclude_fragments: Option<Vec<String>>,

    #[arg(short = 'E', long)]
    /// Use regexp matching instead of literal for keep- and exclude-fragments
    regexp: bool,
}

fn main() -> Result<()> {
    simple_logger::SimpleLogger::new()
        .with_level(LevelFilter::Info)
        .with_colors(true)
        .init()
        .context("failed to initialize simple_logger")?;

    let args = Args::parse();

    let mut results = if !args.inputs.is_empty() {
        RunResult::from_files(&args.inputs)?
    } else {
        log::warn!("Reading results from STDIN");
        RunResult::from_stdin()?
    };

    let mut out: Box<dyn std::io::Write> = if let Some(output) = args.output.as_ref() {
        let out_filename = asgart::utils::make_out_filename(Some(output), "out", &args.format);
        Box::new(File::create(out_filename)?)
    } else {
        Box::new(std::io::stdout())
    };

    let exporter = match args.format.as_str() {
        "json" => Box::new(exporters::JSONExporter) as Box<dyn Exporter>,
        "gff2" => Box::new(exporters::GFF2Exporter) as Box<dyn Exporter>,
        "gff3" => Box::new(exporters::GFF3Exporter) as Box<dyn Exporter>,
        format @ _ => {
            log::warn!("Unknown output format `{}`: using json instead", format);
            Box::new(exporters::JSONExporter) as Box<dyn Exporter>
        }
    };

    if args.collapse {
        results.flatten();
    }
    if args.no_direct {
        results.remove_direct();
    }
    if args.no_reversed {
        results.remove_reversed();
    }
    if args.no_uncomplemented {
        results.remove_uncomplemented();
    }
    if args.no_complemented {
        results.remove_complemented();
    }
    if args.no_inter {
        results.remove_inter();
    }
    if args.no_inter_relaxed {
        results.remove_inter_relaxed();
    }
    if args.no_intra {
        results.remove_intra();
    }
    if let Some(min_length) = args.min_length.as_ref() {
        results.families.iter_mut().for_each(|family| {
            family.retain(|sd| std::cmp::min(sd.left_length, sd.right_length) >= *min_length)
        });
        results.families.retain(|f| !f.is_empty());
    }
    if let Some(max_family_members) = args.max_family_members {
        results.max_family_members(max_family_members);
    }
    if let Some(keep_fragments) = args.keep_fragments.as_ref() {
        if args.regexp {
            for keep_fragment in keep_fragments {
                results
                    .keep_fragments_regexp::<&str>(keep_fragment)
                    .with_context(|| format!("Error while compiling `{}`", keep_fragment))?;
            }
        } else {
            results.keep_fragments(keep_fragments);
        }
    }
    if let Some(restrict_fragments) = args.restrict_fragments.as_ref() {
        if args.regexp {
            for restrict_fragment in restrict_fragments {
                results
                    .restrict_fragments_regexp::<&str>(restrict_fragment)
                    .with_context(|| format!("Error while compiling `{}`", restrict_fragment))?;
            }
        } else {
            results.restrict_fragments(restrict_fragments);
        }
    };
    if let Some(exclude_fragments) = args.exclude_fragments.as_ref() {
        if args.regexp {
            for restrict_fragment in exclude_fragments {
                results
                    .exclude_fragments_regexp::<&str>(restrict_fragment)
                    .with_context(|| format!("Error while compiling `{}`", restrict_fragment))?;
            }
        } else {
            results.exclude_fragments(exclude_fragments);
        }
    }

    exporter.save(&results, &mut out)?;
    Ok(())
}
