use serde_json;
use std::fs::File;
use structs::*;
use std::io::Write;
use ::errors::*;

fn make_filename(basename: &str, ext: &str) -> String {
    if basename.to_lowercase().ends_with(&format!(".{}", ext)) {
        basename.to_owned()
    } else {
        format!("{}.{}", basename, ext)
    }
}

pub trait Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String>;
}

pub struct JSONExporter;
impl Exporter for JSONExporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String> {
        let file_name = make_filename(file_name, "json");
        let mut out = File::create(&file_name).chain_err(|| format!("Unable to create `{}`", &file_name))?;
        let _ = writeln!(&mut out,
                         "{}",
                         serde_json::to_string_pretty(&result).chain_err(|| "Unable to serialize result into JSON")?
        ).chain_err(|| "Unable to write results");

        Ok(file_name)
    }
}


pub struct GFF2Exporter;
impl Exporter for GFF2Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String> {
        let file_name = make_filename(file_name, "gff2");
        let mut out = File::create(&file_name).chain_err(|| format!("Unable to create `{}`", &file_name))?;
        writeln!(&mut out, "track name=Duplications\tuseScore=1\tdescription=\"ASGART - {chr_left} & {chr_right}\"",
                 chr_left = result.strand1.name,
                 chr_right = result.strand2.name
        ).chain_err(|| "Unable to write results")?;
        for (i, family) in result.families.iter().enumerate() {
            for (j, sd) in family.iter().enumerate() {
                writeln!(&mut out,
                         "{chr_left}\tASGART\tSD\t{left}\t{end}\t#{identity}\t#{reverse}\t.\tsd#{i}/{j}-{chr_left}",
                         chr_left = str::replace(&sd.chr_left.trim(), " ", "_"),
                         left     = sd.chr_left_position,
                         end      = sd.chr_left_position + sd.length,
                         identity = sd.identity * 100.0,
                         reverse  = "+",
                         i        = i, j = j
                ).chain_err(|| "Unable to write results")?;
                writeln!(&mut out,
                         "{chr_right}\tASGART\tSD\t{right}\t{end}\t#{identity}\t#{reverse}\t.\tsd#{i}/{j}-{chr_right}",
                         chr_right = str::replace(&sd.chr_right.trim(), " ", "_"),
                         right     = sd.chr_right_position,
                         end       = sd.chr_right_position + sd.length,
                         identity  = sd.identity * 100.0,
                         reverse   = if sd.reversed { "-" } else { "+" },
                         i         = i, j = j
                ).chain_err(|| "Unable to write results")?;
            }
        }

        Ok(file_name)
    }
}



pub struct GFF3Exporter;
impl Exporter for GFF3Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String> {
        let file_name = make_filename(file_name, "gff3");
        let mut out = File::create(&file_name).chain_err(|| format!("Unable to create `{}`", &file_name))?;
        writeln!(&mut out, "##gff-version 3.2.1").chain_err(|| "Unable to write results")?;
        for chr in &result.strand1.map {
            writeln!(
                &mut out,
                "##sequence-region {name} {start} {end}",
                name = &chr.name,
                start  = chr.position + 1,
                end = chr.position + chr.length + 1,
            ).chain_err(|| "Unable to write results")?;
        }
        if result.strand2.name != result.strand1.name {
            for chr in &result.strand2.map {
                writeln!(
                    &mut out,
                    "##sequence-region {name} {start} {end}",
                    name = &chr.name,
                    start  = chr.position + 1,
                    end = chr.position + chr.length + 1,
                ).chain_err(|| "Unable to write results")?;
            }
        }

        for (i, family) in result.families.iter().enumerate() {
            for (j, sd) in family.iter().enumerate() {
                writeln!(
                    &mut out,
                    "{name}\tASGART\tSD\t{left}\t{end}\t{identity}\t{reverse}\t.\tID=sd{i}/{j};Name=SD#{i}/{j}",
                    name     = str::replace(&sd.chr_left.trim(), " ", "_"),
                    left     = sd.chr_left_position + 1,
                    end      = sd.chr_left_position + sd.length + 1,
                    identity = sd.identity,
                    reverse  = "+",
                    i        = i, j = j
                ).chain_err(|| "Unable to write results")?;
                writeln!(
                    &mut out,
                    "{name}\tASGART\tSD\t{right}\t{end}\t{identity}\t{reverse}\t.\tID=sd{i}/{j}-right;Parent=sd{i}/{j};Name=SD#{i}/{j}",
                    name     = str::replace(&sd.chr_right.trim(), " ", "_"),
                    right    = sd.chr_right_position + 1,
                    end      = sd.chr_right_position + sd.length + 1,
                    identity = sd.identity,
                    reverse  = if sd.reversed { "-" } else { "+" },
                    i        = i, j = j
                ).chain_err(|| "Unable to write results")?;
            }
        }

        Ok(file_name)
    }
}
