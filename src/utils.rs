use ::bio::io::fasta;
use std::path;
use ::structs::*;
use ::errors::*;


pub fn complement_nucleotide(n: u8) -> u8 {
    match n {
        b'A' => b'T',
        b'T' => b'A',
        b'G' => b'C',
        b'C' => b'G',
        b'N' => b'N',

        b'a' => b't',
        b't' => b'a',
        b'g' => b'c',
        b'c' => b'g',
        b'n' => b'n',

        _    => b'N'
    }
}


pub fn complemented(text: &[u8]) -> Vec<u8> {
    text.iter().map(|x| complement_nucleotide(*x)).collect::<Vec<u8>>()
}

pub fn slugify(x: &str) -> String {
    x .trim()
        .replace(" ", "_")
        .replace(":", "_")
        .replace("|", "_")
}

pub fn read_fasta(filename: &str, skip_masked: bool) -> Result<(Vec<Start>, Vec<u8>)> {
    let mut map = Vec::new();
    let mut r = Vec::new();

    let reader = fasta::Reader::from_file(filename).chain_err(|| format!("Unable to open `{}`", filename))?;
    let mut counter = 0;

    for record in reader.records() {
        let record = record.chain_err(|| format!("Unable to read {:?}: not a FASTA file", path::Path::new(filename).file_name().unwrap()))?;

        let name = record.id().to_owned();
        let mut seq = record.seq().to_vec();
        if !skip_masked {seq = seq.to_ascii_uppercase();}
        for c in &mut seq {
            if ALPHABET_MASKED.contains(c) && skip_masked {
                *c = b'N'
            } else if !(ALPHABET).contains(c) {
                trace!("Undefined base `{}` replaced by `N`", std::char::from_u32(u32::from(*c)).unwrap());
                *c = b'N'
            }
        }

        map.push(Start {
            name: name,
            position: counter,
            length: seq.len(),
        });
        counter += seq.len();
        r.append(&mut seq);
    }


    Ok((map, r))
}


pub fn make_out_filename(filename: Option<&str>, default: &str, extension: &str) -> std::path::PathBuf {
    filename
        .map(|f| {
            let mut path = std::path::PathBuf::from(f);
            if path.is_dir() { path.push(default); }
            path.set_extension(extension);
            path
        })
        .unwrap_or_else(|| {
                 let mut new = std::path::PathBuf::from(default);
                 new.set_extension(extension);
                 new
        })
}
