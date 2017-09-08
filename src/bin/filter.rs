extern crate asgart;
extern crate rustc_serialize;
extern crate rayon;

use rayon::prelude::*;
use std::env;
use rustc_serialize::json;
use std::io::prelude::*;
use std::io;
use std::io::Write;
use std::fs::File;
use asgart::structs::*;


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
    for chr in &result.strand1.map {
        if chr.name.trim() == name {
            return Some(chr.clone())
        }
    }
    None
}

fn find_right_chromosome(result: &RunResult, name: &str) -> Option<Start> {
    for chr in &result.strand2.map {
        if chr.name.trim() == name {
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

    let mut masks_left = read_masks(masks_left_file).unwrap();
    let mut masks_right = read_masks(masks_right_file).unwrap();

    println!("Precomputing masks...");
    for mask_left in &mut masks_left {
        let chr = find_left_chromosome(&result, &mask_left.chr).expect(&format!("Unable to find `{}`", &mask_left.chr));
        let start = chr.position + mask_left.start;
        let end = chr.position + mask_left.end;
        mask_left.start = start;
        mask_left.end = end;
    }

    for mask_right in &mut masks_right {
        let chr = find_right_chromosome(&result, &mask_right.chr).expect(&format!("Unable to find `{}`", &mask_right.chr));
        let start = chr.position + mask_right.start;
        let end = chr.position + mask_right.end;
        mask_right.start = start;
        mask_right.end = end;
    }
    println!("Done.");

    let old_sds = result.sds.clone();
    result.sds = old_sds.into_par_iter().filter(|sd| {
        let mut keep = true;
        for mask in &masks_left {
            if sd.left >= mask.start && sd.left <= mask.end { keep = false; continue }
        }

        for mask in &masks_right {
            if sd.right >= mask.start && sd.right <= mask.end { keep = false; continue }
        }

        keep
    }).collect();

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
