use std::io::Write;

use anyhow::{Context, Result};
use serde_json;

use crate::structs::*;

pub trait Exporter {
    fn save(&self, result: &RunResult, out: &mut dyn Write) -> Result<()>;
}

pub struct JSONExporter;
impl Exporter for JSONExporter {
    fn save(&self, result: &RunResult, out: &mut dyn Write) -> Result<()> {
        let _ = writeln!(
            out,
            "{}",
            serde_json::to_string_pretty(&result)
                .context("Unable to serialize result into JSON")?
        )
        .context("Unable to write results");

        Ok(())
    }
}

pub struct GFF2Exporter;
impl Exporter for GFF2Exporter {
    fn save(&self, result: &RunResult, out: &mut dyn Write) -> Result<()> {
        writeln!(
            out,
            "track name=Duplications\tuseScore=1\tdescription=\"ASGART - {dataset}\"",
            dataset = result.strand.name,
        )
        .context("Unable to write results")?;
        for (i, family) in result.families.iter().enumerate() {
            for (j, sd) in family.iter().enumerate() {
                writeln!(out,
                         "{chr_left}\tASGART\tSD\t{left}\t{end}\t#{identity}\t+\t.\tSD#{i}/{j}-{chr_left}",
                         chr_left = str::replace(&sd.chr_left.trim(), " ", "_"),
                         left     = sd.chr_left_position,
                         end      = sd.chr_left_position + sd.left_length,
                         identity = sd.identity * 100.0,
                         i        = i, j = j
                ).context("Unable to write results")?;
                writeln!(out,
                         "{chr_right}\tASGART\tSD\t{right}\t{end}\t#{identity}\t#{reverse}\t.\tSD#{i}/{j}-{chr_right}",
                         chr_right = str::replace(&sd.chr_right.trim(), " ", "_"),
                         right     = sd.chr_right_position,
                         end       = sd.chr_right_position + sd.right_length,
                         identity  = sd.identity * 100.0,
                         reverse   = if sd.reversed { "-" } else { "+" },
                         i         = i, j = j
                ).context("Unable to write results")?;
            }
            writeln!(out).context("Unable to write results")?;
        }

        Ok(())
    }
}

pub struct GFF3Exporter;
impl Exporter for GFF3Exporter {
    fn save(&self, result: &RunResult, out: &mut dyn Write) -> Result<()> {
        writeln!(out, "##gff-version 3.2.1").context("Unable to write results")?;
        for chr in &result.strand.map {
            writeln!(
                out,
                "##sequence-region {name} {start} {end}",
                name = &chr.name,
                start = chr.position + 1,
                end = chr.position + chr.length + 1,
            )
            .context("Unable to write results")?;
        }

        for (i, family) in result.families.iter().enumerate() {
            for (j, sd) in family.iter().enumerate() {
                writeln!(
                    out,
                    "{chr_left}\tASGART\tSD\t{left}\t{end}\t{identity}\t+\t.\tID=SD#{i}-{j};Name=SD#{i}-{j}",
                    chr_left  = str::replace(&sd.chr_left.trim(), " ", "_"),
                    left      = sd.chr_left_position + 1,
                    end       = sd.chr_left_position + sd.left_length + 1,
                    identity  = sd.identity,
                    i         = i, j = j
                ).context("Unable to write results")?;
                writeln!(
                    out,
                    "{chr_right}\tASGART\tSD\t{right}\t{end}\t{identity}\t{reverse}\t.\tID=SD#{i}-{j}-right;Parent=SD#{i}-{j};Name=SD#{i}-{j}",
                    chr_right = str::replace(&sd.chr_right.trim(), " ", "_"),
                    right     = sd.chr_right_position + 1,
                    end       = sd.chr_right_position + sd.right_length + 1,
                    identity  = sd.identity,
                    reverse   = if sd.reversed { "-" } else { "+" },
                    i         = i, j = j
                ).context("Unable to write results")?;
            }
            writeln!(out).context("Unable to write results")?;
        }

        Ok(())
    }
}
