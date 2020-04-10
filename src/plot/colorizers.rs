use std::collections::HashMap;

use palette::*;
use rand::prelude::*;

use crate::structs::*;

pub trait Colorizer {
    fn color_fragment(&self, name: &str) -> String;
    fn color(&self, d: &SD) -> String;
}

pub struct TypeColorizer {
    direct_color: Srgb,
    rc_color: Srgb,
}
impl TypeColorizer {
    pub fn new(direct_color: (f32, f32, f32), rc_color: (f32, f32, f32)) -> Self {
        TypeColorizer {
            direct_color: Srgb::from_components(direct_color),
            rc_color: Srgb::from_components(rc_color),
        }
    }
}
impl Colorizer for TypeColorizer {
    fn color_fragment(&self, _name: &str) -> String {
        "#cccccc".to_owned()
    }
    fn color(&self, sd: &SD) -> String {
        let color = (if !sd.reversed && !sd.complemented {
            self.direct_color
        } else {
            self.rc_color
        })
        .into_components();

        format!(
            "#{:02x}{:02x}{:02x}",
            (color.0 * 255.0) as u8,
            (color.1 * 255.0) as u8,
            (color.2 * 255.0) as u8
        )
    }
}

pub struct PositionColorizer {
    total_length: f64,
    gradient: Gradient<Hsv<encoding::Srgb, f64>>,
}
impl PositionColorizer {
    pub fn new(result: &RunResult) -> Self {
        PositionColorizer {
            total_length: result.strand.length as f64,
            gradient: Gradient::new(vec![
                Hsv::from(LinSrgb::new(1.0, 0.1, 0.1)),
                Hsv::from(LinSrgb::new(0.1, 1.0, 1.0)),
            ]),
        }
    }
}
impl Colorizer for PositionColorizer {
    fn color_fragment(&self, _name: &str) -> String {
        "#cccccc".to_owned()
    }
    fn color(&self, sd: &SD) -> String {
        let color = Srgb::from_hsv(
            self.gradient
                .get(sd.global_left_position as f64 / self.total_length),
        )
        .into_components();

        format!(
            "#{:02x}{:02x}{:02x}",
            (color.0 * 255.0) as u8,
            (color.1 * 255.0) as u8,
            (color.2 * 255.0) as u8
        )
    }
}

pub struct FragmentColorizer {
    colors: HashMap<String, rgb::Rgb<encoding::Srgb, f64>>,
}
impl FragmentColorizer {
    pub fn new(result: &RunResult) -> Self {
        let mut colors = (0..result.strand.map.len())
            .map(|i| {
                rgb::Rgb::from_hsv(Hsv::new(
                    30.0 + 330.0 * (i as f64) / (result.strand.map.len() as f64),
                    1.0,
                    0.7,
                ))
            })
            .collect::<Vec<_>>();
        colors.shuffle(&mut thread_rng());

        let names = result
            .strand
            .map
            .iter()
            .map(|chr| chr.name.clone())
            .collect::<Vec<String>>();

        FragmentColorizer {
            colors: names.into_iter().zip(colors.into_iter()).collect(),
        }
    }
}
impl Colorizer for FragmentColorizer {
    fn color_fragment(&self, name: &str) -> String {
        let color = self
            .colors
            .get(name)
            .expect(&format!("Unable to get {}", name))
            .into_components();
        format!(
            "#{:02x}{:02x}{:02x}",
            (color.0 / 1.3 * 255.0) as u8,
            (color.1 / 1.3 * 255.0) as u8,
            (color.2 / 1.3 * 255.0) as u8
        )
    }

    fn color(&self, sd: &SD) -> String {
        let color = self
            .colors
            .get(&sd.chr_left)
            .expect(&format!("Unable to get {}", &sd.chr_left))
            .into_components();

        format!(
            "#{:02x}{:02x}{:02x}",
            (color.0 * 255.0) as u8,
            (color.1 * 255.0) as u8,
            (color.2 * 255.0) as u8
        )
    }
}
