extern crate palette;
extern crate separator;

use ::structs::*;

pub mod chord_plot;
pub mod flat_plot;

pub struct Settings {
    pub result: RunResult,
    pub out_file: String,

    pub size: f64,

    pub thickness: f64,
    pub color1: String,
    pub color2: String,

    pub min_length: usize,
}

pub trait Plotter {
    fn new(settings: Settings) -> Self;
    fn plot(self);
}
