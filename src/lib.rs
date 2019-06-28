#[macro_use] extern crate serde_derive;
extern crate superslice;
extern crate serde;
extern crate serde_json;
extern crate separator;
extern crate rayon;
#[macro_use] extern crate error_chain;

pub mod structs;
pub mod logger;
pub mod utils;
pub mod divsufsort64;
pub mod searcher;
pub mod automaton;
pub mod plot;
pub mod exporters;
pub mod errors {
    error_chain!{}
}
