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
        writeln!(&mut out, "track name=Duplications\tuseScore=1\tdescription=\"ASGART - {chr_left} vs. {chr_right}\"",
                 chr_left = result.strand1.name,
                 chr_right = result.strand2.name
        ).chain_err(|| "Unable to write results")?;
        for (i, sd) in result.sds.iter().enumerate() {
            writeln!(&mut out,
                     "{chr_left}\tASGART\tSD\t{left}\t{end}\t#{identity}\t#{reverse}\t.\tsd#{i}-{chr_left}",
                     chr_left = result.strand1.name,
                     left = sd.left,
                     end = sd.left + sd.length,
                     identity = sd.identity * 10.0,
                     reverse = "+",
                     i = i
            ).chain_err(|| "Unable to write results")?;
            writeln!(&mut out,
                     "{chr_right}\tASGART\tSD\t{right}\t{end}\t#{identity}\t#{reverse}\t.\tsd#{i}-{chr_right}",
                     chr_right = result.strand2.name,
                     right = sd.right,
                     end = sd.right + sd.length,
                     identity = sd.identity * 10.0,
                     reverse = if sd.reversed { "-" } else { "+" },
                     i = i
            ).chain_err(|| "Unable to write results")?;
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
            let start = chr.position + 1;
            let end = chr.position + chr.length + 1;
            let name = &chr.name;
            writeln!(&mut out, "##sequence-region {} {} {}", name, start, end).chain_err(|| "Unable to write results")?;
        }
        for chr in &result.strand2.map {
            let start = chr.position + 1;
            let end = chr.position + chr.length + 1;
            let name = &chr.name;
            writeln!(&mut out, "##sequence-region {} {} {}", name, start, end).chain_err(|| "Unable to write results")?;
        }

        for (i, sd) in result.sds.iter().enumerate() {
            let left_chr = result.strand1.find_chr_by_pos(sd.left+1);
            let right_chr = result.strand2.find_chr_by_pos(sd.right+1);

            writeln!(&mut out,
                     "{name}\tASGART\tSD\t{left}\t{end}\t{identity}\t{reverse}\t.\tID=sd{i};Name=SD#{i}",
                     name     = str::replace(&left_chr.name.trim(), " ", "_"),
                     left     = sd.left + 1 - left_chr.position,
                     end      = sd.left + 1 + sd.length - left_chr.position,
                     identity = sd.identity,
                     reverse  = "+",
                     i        = i
            ).chain_err(|| "Unable to write results")?;
            writeln!(&mut out,
                     "{name}\tASGART\tSD\t{right}\t{end}\t{identity}\t{reverse}\t.\tID=sd{i}-right;Parent=sd{i};Name=SD#{i}",
                     name     = str::replace(&right_chr.name.trim(), " ", "_"),
                     right    = sd.right + 1 - right_chr.position,
                     end      = sd.right + sd.length + 1 - right_chr.position,
                     identity = sd.identity,
                     reverse  = if sd.reversed { "-" } else { "+" },
                     i        = i
            ).chain_err(|| "Unable to write results")?;
        }

        Ok(file_name)
    }
}
