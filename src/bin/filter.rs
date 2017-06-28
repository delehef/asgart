extern crate bio;
extern crate asgart;
extern crate rustc_serialize;

use std::process::Command;
use asgart::utils;
use std::path;
use std::env;
use rustc_serialize::json;
use std::io::prelude::*;
use std::io;
use std::io::Write;
use std::fs::File;
use std::ascii::AsciiExt;
use asgart::structs::*;
use bio::io::fasta;

fn old_filter(json_file: &str) {
    let mut f = File::open(json_file).unwrap();
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    let result: RunResult = json::decode(&s).unwrap();

    let strand1 = read_fasta(&result.strand1.name).unwrap();
    let strand2 = if result.strand1.name == result.strand2.name {
        strand1.clone()
    } else {
        read_fasta(&result.strand2.name).unwrap()
    };

    for sd in result.sds {
        let mut fasta_right = strand2[sd.left .. sd.left + sd.length].to_vec();
        if sd.reversed { fasta_right.reverse(); }
        if sd.translated { fasta_right = utils::translated(&fasta_right); }

        File::create("sd_left.fasta")
            .expect("Unable to create sd_left.fasta")
            .write_all(&strand1[sd.left .. sd.left + sd.length]);
        File::create("sd_right.fasta")
            .expect("Unable to create sd_right.fasta")
            .write_all(&fasta_right);

        println!("/home/franklin/repos/reputer/repfind  -p -l 1000 sd_left.fasta sd_right.fasta");
        let output = Command::new("/home/franklin/repos/reputer/repfind")
            .arg("-p")
            .arg("-l 10")
            .arg("sd_left.fasta")
            .arg("sd_right.fasta")
            .output()
            .expect("Failed to launch REPuter")
            .stdout;

        println!("{}", std::str::from_utf8(&output).unwrap());
    }
}

fn read_fasta(filename: &str) -> Result<Vec<u8>, io::Error> {
    let mut r = Vec::new();
    let reader = fasta::Reader::from_file(filename).expect("Unable to read");

    for record in reader.records() {
        let record = record.expect(&format!("Unable to read {:?}: not a FASTA file", path::Path::new(filename).file_name().unwrap()));
        r.append(&mut record.seq().to_vec().to_ascii_uppercase());
    }

    Ok(r)
}

struct Mask {
    chr: String,
    start: usize,
    end: usize,
}

fn read_masks(masks_file: &str) -> Result<Vec<Mask>, io::Error> {
    let mut r = Vec::new();
    let mut f = try!(File::open(masks_file));
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    let masks_string = str::replace(&s, ">", "");

    for l in masks_string.lines() {
        let v : Vec<&str> = l.split('\t').collect();
        assert_eq!(v.len(), 3, "Error");
        r.push(Mask {
            chr:   v[0].to_owned(),
            start: v[1].parse::<usize>().unwrap(),
            end:   v[2].parse::<usize>().unwrap()
        });
    }

    Ok(r)
}

fn find_left_chromosome(result: &RunResult, name: &str) -> Option<Start> {
    for chr in result.strand1.map.iter() {
        if chr.name == name {
            return Some(chr.clone())
        }
    }
    None
}

fn find_right_chromosome(result: &RunResult, name: &str) -> Option<Start> {
    for chr in result.strand2.map.iter() {
        if chr.name == name {
            return Some(chr.clone())
        }
    }
    None
}

fn filter(json_file: &str, masks_left_file: &str, masks_right_file: &str) -> RunResult {
    let mut f = File::open(json_file).expect(&format!("Unable to open {}", json_file));
    let mut s = String::new();
    let _ = f.read_to_string(&mut s);
    let mut result: RunResult = json::decode(&s).unwrap();
    let mut new_sds : Vec<SD> = Vec::new();

    let mut masks_left = read_masks(masks_left_file).unwrap();
    let mut masks_right = read_masks(masks_right_file).unwrap();

    println!("Precomputing masks...");
    for mask_left in masks_left.iter_mut() {
        let chr = find_left_chromosome(&result, &mask_left.chr).expect(&format!("Unable to find `{}`", &mask_left.chr));
        let start = chr.position + mask_left.start;
        let end = chr.position + mask_left.end;
        mask_left.start = start;
        mask_left.end = end;
    }

    for mask_right in masks_right.iter_mut() {
        let chr = find_right_chromosome(&result, &mask_right.chr).expect(&format!("Unable to find `{}`", &mask_right.chr));
        let start = chr.position + mask_right.start;
        let end = chr.position + mask_right.end;
        mask_right.start = start;
        mask_right.end = end;
    }
    println!("Done.");

    let mut i = 0;
    for sd in result.sds.iter() {
        i += 1;
        println!("{}/{}", i, result.sds.len());
        let mut keep = true;
        for mask in masks_left.iter() {
            if sd.left >= mask.start && sd.left <= mask.end { keep = false; continue }
        }

        for mask in masks_right.iter() {
            if sd.right >= mask.start && sd.right <= mask.end { keep = false; continue }
        }

        if keep { println!("Keeping"); new_sds.push(sd.clone()) }
    }

   result.sds = new_sds; 

   result
}

fn main() {
    let args : Vec<String> = env::args().skip(1).collect();
    let json_file = &args[0];
    let masks_left_file = &args[2];
    let masks_right_file = &args[2];
    println!("Working on {} & {}, {}", json_file, masks_left_file, masks_right_file);
    let result = filter(json_file, masks_left_file, masks_right_file);
    let mut out = File::create(&format!("{}.masked", json_file)).expect("Unable to create out file");
    println!("Final: {}", result.sds.len());
    writeln!(&mut out,
             "{}",
             rustc_serialize::json::as_pretty_json(&result)).expect("Unable to write results");
}
