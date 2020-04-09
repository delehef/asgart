pub extern crate log;
#[macro_use]
pub extern crate serde_derive;
pub extern crate bio;
pub extern crate colored;
pub extern crate console;
pub extern crate indicatif;
pub extern crate num_cpus;
pub extern crate rand;
pub extern crate rayon;
pub extern crate regex;
pub extern crate serde;
pub extern crate serde_json;
pub extern crate superslice;
pub extern crate threadpool;
pub extern crate thousands;
#[macro_use]
pub extern crate error_chain;

pub mod structs;
#[macro_use]
pub mod logger;
pub mod automaton;
pub mod divsufsort;
pub mod exporters;
pub mod plot;
pub mod searcher;
pub mod utils;
pub mod errors {
    error_chain! {}
}
