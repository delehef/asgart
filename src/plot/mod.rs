use anyhow::{Context, Result};

use crate::structs::*;

pub mod chord_plot;
pub mod flat_plot;
// pub mod eye_plot;
pub mod circos_plot;
pub mod colorizers;
pub mod genome_plot;
pub mod rosary_plot;

pub struct Settings {
    pub out_file: String,
    pub size: f64,

    pub min_thickness: f64,
    pub color1: String,
    pub color2: String,

    pub feature_tracks: Vec<Vec<Feature>>,
}

#[derive(Debug, Clone)]
pub enum FeaturePosition {
    Relative {
        chr: String,
        start: usize,
        length: usize,
    },
    Absolute {
        start: usize,
        length: usize,
    },
}

#[derive(Debug, Clone)]
pub struct Feature {
    pub name: String,
    pub positions: Vec<FeaturePosition>,
}

pub trait Plotter {
    fn plot(&self) -> Result<()>;
}

pub enum SvgObject {
    Line {
        x1: f32,
        y1: f32,
        x2: f32,
        y2: f32,
        stroke: Option<String>,
        stroke_width: f32,
        hover: Option<String>,
    },

    Circle {
        cx: f32,
        cy: f32,
        r: f32,
        fill: String,
    },

    Text {
        x: f32,
        y: f32,
        text: String,
        font_size: Option<f32>,
        color: Option<String>,
    },
}

#[derive(Debug, Clone, Copy)]
pub struct BBox {
    x1: f32,
    y1: f32,
    x2: f32,
    y2: f32,
}

impl SvgObject {
    fn render(&self) -> String {
        match self {
            SvgObject::Line {
                x1,
                y1,
                x2,
                y2,
                stroke,
                stroke_width,
                hover,
            } => {
                let mut style = format!("stroke-width: {};", stroke_width);
                if let Some(stroke) = stroke {
                    style.push_str(&format!("stroke: {};", stroke));
                }
                let mut inner = format!("x1='{}' y1='{}' x2='{}' y2='{}'", x1, y1, x2, y2);
                if !style.is_empty() {
                    inner.push_str(&format!(" style='{}'", style));
                }
                match hover {
                    Some(text) => format!("<line {}><title>{}</title></line>", inner, text),
                    None => format!("<line {}/>", inner),
                }
            }
            SvgObject::Circle { cx, cy, r, fill } => format!(
                "<circle cx='{}' cy='{}' r='{}' fill='{}'/>",
                cx, cy, r, fill
            ),
            SvgObject::Text {
                x,
                y,
                text,
                font_size,
                color,
            } => format!(
                "<text x='{}' y='{}' font-family='Helvetica' fill='{}' font-size='{}'>{}</text>",
                x,
                y,
                color.as_ref().unwrap_or(&"#000".to_string()),
                font_size.unwrap_or(10.0),
                text
            ),
        }
    }

    pub fn shift(&mut self, dx: f32, dy: f32) {
        match self {
            SvgObject::Line {
                ref mut x1,
                ref mut y1,
                ref mut x2,
                ref mut y2,
                ..
            } => {
                *x1 += dx;
                *x2 += dx;
                *y1 += dy;
                *y2 += dy;
            }
            SvgObject::Circle {
                ref mut cx,
                ref mut cy,
                ..
            } => {
                *cx += dx;
                *cy += dy;
            }
            SvgObject::Text {
                ref mut x,
                ref mut y,
                ..
            } => {
                *x += dx;
                *y += dy;
            }
        };
    }

    pub fn scale(&mut self, s: f32) {
        match *self {
            SvgObject::Line {
                ref mut x1,
                ref mut y1,
                ref mut x2,
                ref mut y2,
                ref mut stroke_width,
                ..
            } => {
                *x1 *= s;
                *x2 *= s;
                *y1 *= s;
                *y2 *= s;
                *stroke_width *= s;
            }
            SvgObject::Circle {
                ref mut cx,
                ref mut cy,
                ref mut r,
                ..
            } => {
                *cx *= s;
                *cy *= s;
                *r *= s;
            }
            SvgObject::Text {
                ref mut x,
                ref mut y,
                ref mut font_size,
                ..
            } => {
                *x *= s;
                *y *= s;
                font_size.map(|x| x * s);
            }
        }
    }

    pub fn dims(&self) -> (f32, f32) {
        match self {
            SvgObject::Line { x1, y1, x2, y2, .. } => ((x2 - x1).abs(), (y2 - y1).abs()),
            SvgObject::Circle { r, .. } => (2. * r, 2. * r),
            SvgObject::Text {
                font_size, text, ..
            } => (
                font_size.unwrap_or(10.) * text.len() as f32,
                font_size.unwrap_or(10.),
            ),
        }
    }

    pub fn bbox(&self) -> BBox {
        match self {
            SvgObject::Line {
                x1,
                y1,
                x2,
                y2,
                stroke_width,
                ..
            } => {
                let (x_min, x_max) = if x1 > x2 { (x2, x1) } else { (x1, x2) };
                let (y_min, y_max) = if y1 > y2 { (y2, y1) } else { (y1, y2) };

                // FIXME ugly approximation
                BBox {
                    x1: x_min - stroke_width / 2.,
                    y1: y_min - stroke_width / 2.,
                    x2: x_max + stroke_width / 2.,
                    y2: y_max + stroke_width / 2.,
                }
            }
            SvgObject::Circle { cx, cy, r, .. } => BBox {
                x1: cx - r,
                y1: cy - r,
                x2: cx + r,
                y2: cy + r,
            },
            SvgObject::Text {
                x,
                y,
                font_size,
                text,
                ..
            } => BBox {
                x1: *x,
                y1: *y,
                x2: x + font_size.unwrap_or(10.) * text.len() as f32,
                y2: y + font_size.unwrap_or(10.),
            },
        }
    }

    pub fn transpose(&mut self) {
        match self {
            SvgObject::Line {
                ref mut x1,
                ref mut y1,
                ref mut x2,
                ref mut y2,
                ..
            } => {
                std::mem::swap(x1, y1);
                std::mem::swap(x2, y2);
            }
            SvgObject::Circle {
                ref mut cx,
                ref mut cy,
                ..
            } => {
                std::mem::swap(cx, cy);
            }
            SvgObject::Text {
                ref mut x,
                ref mut y,
                ..
            } => {
                std::mem::swap(x, y);
            }
        }
    }
}

pub struct SvgGroup {
    content: Vec<SvgObject>,
}
impl Default for SvgGroup {
    fn default() -> Self {
        Self::new()
    }
}

impl SvgGroup {
    pub fn new() -> Self {
        SvgGroup {
            content: Vec::new(),
        }
    }

    pub fn push(mut self, o: SvgObject) -> Self {
        self.content.push(o);
        self
    }

    pub fn append(mut self, other: Self) -> Self {
        self.content.extend(other.content);
        self
    }

    pub fn extend(mut self, other: impl Iterator<Item = SvgObject>) -> Self {
        self.content.extend(other);
        self
    }

    pub fn render(&self) -> String {
        self.content
            .iter()
            .map(SvgObject::render)
            .collect::<Vec<String>>()
            .join("\n")
    }

    pub fn shift(mut self, dx: f32, dy: f32) -> Self {
        self.content.iter_mut().for_each(|o| o.shift(dx, dy));
        self
    }

    pub fn scale(mut self, s: f32) -> Self {
        self.content.iter_mut().for_each(|o| o.scale(s));
        self
    }

    pub fn bbox(&self) -> BBox {
        let (mut x1, mut y1, mut x2, mut y2) = (0., 0., 0., 0.);
        for x in self.content.iter() {
            let bbox = x.bbox();

            if bbox.x1 < x1 {
                x1 = bbox.x1
            };
            if bbox.y1 < y1 {
                y1 = bbox.y1
            };
            if bbox.x2 > x2 {
                x2 = bbox.x2
            };
            if bbox.y2 > y2 {
                y2 = bbox.y2
            };
        }

        BBox { x1, y1, x2, y2 }
    }

    pub fn dims(&self) -> (f32, f32) {
        let bbox = self.bbox();
        (bbox.x2 - bbox.x1, bbox.y2 - bbox.y1)
    }

    pub fn transpose(mut self) -> Self {
        self.content.iter_mut().for_each(|x| x.transpose());
        self
    }
}
