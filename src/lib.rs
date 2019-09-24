#[macro_use] extern crate serde_derive;
extern crate superslice;
extern crate serde;
pub extern crate serde_json;
pub extern crate separator;
pub extern crate rayon;
#[macro_use] extern crate error_chain;
#[macro_use] pub extern crate log;

pub mod structs;
pub mod logger;
pub mod utils;
pub mod divsufsort;
pub mod searcher;
pub mod automaton;
pub mod plot;
pub mod exporters;
pub mod errors {
    error_chain!{}
}
