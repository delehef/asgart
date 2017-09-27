extern crate rustc_serialize;
extern crate separator;

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
