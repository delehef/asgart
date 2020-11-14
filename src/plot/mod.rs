use anyhow::{Context, Result};

use crate::structs::*;

pub mod chord_plot;
pub mod flat_plot;
// pub mod eye_plot;
pub mod circos_plot;
pub mod colorizers;
pub mod genome_plot;

pub struct Settings {
    pub out_file: String,
    pub size:     f64,

    pub min_thickness: f64,
    pub color1:        String,
    pub color2:        String,

    pub feature_tracks: Vec<Vec<Feature>>,
}

#[derive(Debug, Clone)]
pub enum FeaturePosition {
    Relative {
        chr:    String,
        start:  usize,
        length: usize,
    },
    Absolute {
        start:  usize,
        length: usize,
    },
}

#[derive(Debug, Clone)]
pub struct Feature {
    pub name:      String,
    pub positions: Vec<FeaturePosition>,
}

pub trait Plotter {
    fn plot(&self) -> Result<()>;
}


pub enum SvgObject {
    Line {
        x1: f32, y1: f32,
        x2: f32, y2: f32,
        stroke: Option<String>,
        stroke_width: f32,
        hover: Option<String>,
    },

    Circle {
        cx: f32, cy: f32, r: f32,
        fill: String,
    }
}

impl SvgObject {
    fn render(&self) -> String {
        match self {
            SvgObject::Line { x1, y1, x2, y2, stroke, stroke_width, hover } => {
                let mut style = format!("stroke-width: {};", stroke_width);
                if let Some(stroke) = stroke {
                    style.push_str(&format!("stroke: {};", stroke));
                }
                let mut inner = format!("x1='{}' y1='{}' x2='{}' y2='{}'", x1, y1, x2, y2);
                if !style.is_empty() {
                    inner.push_str(&format!(" style='{}'", style));
                }
                match hover {
                    Some(text) => {
                        format!("<line {}><title>{}</title></line>", inner, text)
                    }
                    None => {
                        format!("<line {}/>", inner)
                    }
                }
            },
            SvgObject::Circle { cx, cy, r, fill } => {
                format!("<circle cx='{}' cy='{}' r='{}' fill='{}'/>", cx, cy, r, fill)
            },
        }
    }
}

pub struct SvgGroup {
    content: Vec<SvgObject>,
}
impl SvgGroup {
    pub fn new() -> Self {
        SvgGroup { content: Vec::new() }
    }

    pub fn push(&mut self, o: SvgObject) {
        self.content.push(o)
    }

    pub fn render(&self) -> String {
        self.content
            .iter()
            .map(SvgObject::render)
            .collect::<Vec<String>>()
            .join("\n")
    }
}
