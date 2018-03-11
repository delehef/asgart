use ::structs::*;

pub mod chord_plot;
pub mod flat_plot;
pub mod eye_plot;

pub struct Settings {
    pub result: RunResult,
    pub out_file: String,

    pub size: f64,

    pub thickness: f64,
    pub color1: String,
    pub color2: String,

    pub min_length: usize,
    pub min_identity: f32,
    pub plot_if_reversed: bool,
    pub plot_if_translated: bool,

    pub gene_tracks: Vec<Vec<Gene>>,
}

#[derive(Debug, Clone)]
pub enum GenePosition {
    Relative { chr: String, start: usize, length: usize },
    Absolute { start: usize, length: usize }
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
