use rustc_serialize;
use std::fs::File;
use structs::*;
// use asgart::errors::*;
use std::io::Write;
use ::errors::*;


pub trait Exporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String>;
}

pub struct JSONExporter;
impl Exporter for JSONExporter {
    fn save(&self, result: &RunResult, file_name: &str) -> Result<String> {
        let file_name = if file_name.to_lowercase().ends_with(".json") {
            file_name.to_owned()
        } else {
            format!("{}.json", file_name)
        };
        let mut out = File::create(&file_name).chain_err(|| format!("Unable to create `{}`", &file_name))?;
        writeln!(&mut out,
                 "{}",
                 rustc_serialize::json::as_pretty_json(&result)).expect("Unable to write results");

        Ok(file_name)
    }
}
