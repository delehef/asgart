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

    pub gene_tracks: Vec<Gene>,
}

#[derive(Debug, Clone)]
pub enum GenePosition {
    Relative { chr: String, position: usize },
    Absolute { position: usize }
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub name: String,
    pub positions: Vec<GenePosition>,
}

pub trait Plotter {
    fn new(settings: Settings) -> Self;
    fn plot(self);
}
