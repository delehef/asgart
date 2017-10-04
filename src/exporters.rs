use rustc_serialize;
use std::fs::File;
use structs::*;
// use asgart::errors::*;
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
        writeln!(&mut out, "{}", rustc_serialize::json::as_pretty_json(&result)).expect("Unable to write results");

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
        ).expect("Unable to write results");
        for (i, sd) in result.sds.iter().enumerate() {
            writeln!(&mut out,
                     "{chr_left}\tASGART\tSD\t{left}\t{end}\t#{identity}\t#{reverse}\t.\tsd#{i}-{chr_left}",
                     chr_left = result.strand1.name,
                     left = sd.left,
                     end = sd.left + sd.length,
                     identity = sd.identity * 10.0,
                     reverse = if sd.reversed { "-" } else { "+" },
                     i = i
            ).expect("Unable to write results");
            writeln!(&mut out,
                     "{chr_right}\tASGART\tSD\t{right}\t{end}\t#{identity}\t#{reverse}\t.\tsd#{i}-{chr_right}",
                     chr_right = result.strand2.name,
                     right = sd.right,
                     end = sd.right + sd.length,
                     identity = sd.identity * 10.0,
                     reverse = if sd.reversed { "-" } else { "+" },
                     i = i
            ).expect("Unable to write results");
        }

        Ok(file_name)
    }
}



pub struct GFF3Exporter;
impl Exporter for GFF3Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String> {
        let file_name = make_filename(file_name, "gff3");
        let mut out = File::create(&file_name).chain_err(|| format!("Unable to create `{}`", &file_name))?;
        writeln!(&mut out, "##gff-version 3.2.1").expect("Unable to write results");
        for chr in &result.strand1.map {
            let start = chr.position + 1;
            let end = chr.position + chr.length + 1;
            let name = &chr.name;
            writeln!(&mut out, "##sequence-region {} {} {}", name, start, end).expect("Unable to write results");
        }
        for chr in &result.strand2.map {
            let start = chr.position + 1;
            let end = chr.position + chr.length + 1;
            let name = &chr.name;
            writeln!(&mut out, "##sequence-region {} {} {}", name, start, end).expect("Unable to write results");
        }

        for (i, sd) in result.sds.iter().enumerate() {
            writeln!(&mut out,
                     "{name}\tASGART\t.\t{left}\t{end}\t{identity}\t+\t.\tID=sd{i}",
                     name = str::replace(&result.strand1.map[0].name, " ", "_"),
                     left = sd.left+1,
                     end = sd.left + 1 + sd.length,
                     identity = sd.identity,
                     i = i
            ).expect("Unable to write results");
            writeln!(&mut out,
                     "{name}\tASGART\t.\t{right}\t{end}\t{identity}\t+\t.\tID=sd{i}-right;Parent=sd{i}",
                     name = str::replace(&result.strand2.map[0].name, " ", "_"),
                     right = sd.right + 1,
                     end = sd.right + sd.length + 1,
                     identity = sd.identity,
                     i = i
            ).expect("Unable to write results");
        }

        Ok(file_name)
    }
}
