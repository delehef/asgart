#[macro_use] pub extern crate log;
#[macro_use] pub extern crate serde_derive;
pub extern crate superslice;
pub extern crate serde;
pub extern crate serde_json;
pub extern crate separator;
pub extern crate rayon;
pub extern crate bio;
pub extern crate colored;
pub extern crate regex;
pub extern crate rand;
pub extern crate threadpool;
pub extern crate indicatif;
pub extern crate console;
pub extern crate num_cpus;
#[macro_use] pub extern crate error_chain;

pub mod structs;
#[macro_use] pub mod logger;
pub mod utils;
pub mod divsufsort;
pub mod searcher;
pub mod automaton;
pub mod plot;
pub mod exporters;
pub mod errors {
    error_chain!{}
}
