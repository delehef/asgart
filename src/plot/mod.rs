extern crate regex;

use ::errors::*;
use ::structs::*;

pub mod chord_plot;
pub mod flat_plot;
pub mod eye_plot;
pub mod genome_plot;
pub mod circos_plot;

pub struct Settings {
    pub out_file: String,

    pub size: f64,

    pub thickness: f64,
    pub color1: String,
    pub color2: String,

    pub min_length: usize,
    pub min_identity: f32,
    pub filter_direct: bool,
    pub filter_reversed: bool,
    pub filter_complemented: bool,
    pub filter_non_complemented: bool,

    pub feature_tracks: Vec<Vec<Feature>>,
}

#[derive(Debug, Clone)]
pub enum FeaturePosition {
    Relative { chr: String, start: usize, length: usize },
    Absolute { start: usize, length: usize }
}

#[derive(Debug, Clone)]
pub struct Feature {
    pub name: String,
    pub positions: Vec<FeaturePosition>,
}

pub trait Plotter {
    fn new(settings: Settings, result: RunResult) -> Self;
    fn plot(&self) -> Result<()>;
}
