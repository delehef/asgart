use serde_json;
use std::fs::File;
use structs::*;
use std::io::Write;
use ::errors::*;

pub trait Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<()>;
}

pub struct JSONExporter;
impl Exporter for JSONExporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<()> {
        let mut out = File::create(&file_name)
            .chain_err(|| format!("Unable to create `{}`", file_name))?;


        let _ = writeln!(&mut out,
                         "{}",
                         serde_json::to_string_pretty(&result).chain_err(|| "Unable to serialize result into JSON")?
        ).chain_err(|| "Unable to write results");

        Ok(())
    }
}


pub struct GFF2Exporter;
impl Exporter for GFF2Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<()> {
        let mut out = File::create(&file_name)
            .chain_err(|| format!("Unable to create `{}`", &file_name))?;

        writeln!(&mut out, "track name=Duplications\tuseScore=1\tdescription=\"ASGART - {dataset}\"",
                 dataset = result.strand.name,
        ).chain_err(|| "Unable to write results")?;
        for (i, family) in result.families.iter().enumerate() {
            for (j, sd) in family.iter().enumerate() {
                writeln!(&mut out,
                         "{chr_left}\tASGART\tSD\t{left}\t{end}\t#{identity}\t+\t.\tSD#{i}/{j}-{chr_left}",
                         chr_left = str::replace(&sd.chr_left.trim(), " ", "_"),
                         left     = sd.chr_left_position,
                         end      = sd.chr_left_position + sd.left_length,
                         identity = sd.identity * 100.0,
                         i        = i, j = j
                ).chain_err(|| "Unable to write results")?;
                writeln!(&mut out,
                         "{chr_right}\tASGART\tSD\t{right}\t{end}\t#{identity}\t#{reverse}\t.\tSD#{i}/{j}-{chr_right}",
                         chr_right = str::replace(&sd.chr_right.trim(), " ", "_"),
                         right     = sd.chr_right_position,
                         end       = sd.chr_right_position + sd.right_length,
                         identity  = sd.identity * 100.0,
                         reverse   = if sd.reversed { "-" } else { "+" },
                         i         = i, j = j
                ).chain_err(|| "Unable to write results")?;
            }
            writeln!(&mut out).chain_err(|| "Unable to write results")?;
        }

        Ok(())
    }
}



pub struct GFF3Exporter;
impl Exporter for GFF3Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<()> {
        let mut out = File::create(&file_name)
            .chain_err(|| format!("Unable to create `{}`", &file_name))?;

        writeln!(&mut out, "##gff-version 3.2.1").chain_err(|| "Unable to write results")?;
        for chr in &result.strand.map {
            writeln!(
                &mut out,
                "##sequence-region {name} {start} {end}",
                name = &chr.name,
                start  = chr.position + 1,
                end = chr.position + chr.length + 1,
            ).chain_err(|| "Unable to write results")?;
        }

        for (i, family) in result.families.iter().enumerate() {
            for (j, sd) in family.iter().enumerate() {
                writeln!(
                    &mut out,
                    "{chr_left}\tASGART\tSD\t{left}\t{end}\t{identity}\t+\t.\tID=SD#{i}-{j};Name=SD#{i}-{j}",
                    chr_left  = str::replace(&sd.chr_left.trim(), " ", "_"),
                    left      = sd.chr_left_position + 1,
                    end       = sd.chr_left_position + sd.left_length + 1,
                    identity  = sd.identity,
                    i         = i, j = j
                ).chain_err(|| "Unable to write results")?;
                writeln!(
                    &mut out,
                    "{chr_right}\tASGART\tSD\t{right}\t{end}\t{identity}\t{reverse}\t.\tID=SD#{i}-{j}-right;Parent=SD#{i}-{j};Name=SD#{i}-{j}",
                    chr_right = str::replace(&sd.chr_right.trim(), " ", "_"),
                    right     = sd.chr_right_position + 1,
                    end       = sd.chr_right_position + sd.right_length + 1,
                    identity  = sd.identity,
                    reverse   = if sd.reversed { "-" } else { "+" },
                    i         = i, j = j
                ).chain_err(|| "Unable to write results")?;
            }
            writeln!(&mut out).chain_err(|| "Unable to write results")?;
        }

        Ok(())
    }
}
