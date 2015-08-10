#![feature(collections)]
extern crate bio;

// Import some modules
use bio::alphabets;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::bwt::bwt;
use bio::data_structures::fmindex::FMIndex;
use bio::io::fasta;

fn main ()
{
    let reader = fasta::Reader::from_file("/Users/franklin/test.fasta");
    for record in reader.unwrap().records() {
        let seq = record.unwrap();
        println!("{:?}", suffix_array(seq.seq()).len());
    }
}

