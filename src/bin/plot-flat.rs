#![feature(step_by)]

extern crate asgart;
extern crate rustc_serialize;
extern crate palette;

use std::cmp;
use std::env;
use std::io::prelude::*;
use std::fs::File;
use rustc_serialize::json;
use std::f64::consts::PI;
use palette::{Rgb, Hsv, RgbHue, FromColor};

use asgart::structs::*;

struct FlatPlotter {
    result: RunResult,
    max_length: f64,
    width: f64,
    height: f64,
    svg: String,
}

impl FlatPlotter {
    fn new(json_file: &str, width: f64) -> FlatPlotter {
        let mut f = File::open(json_file).unwrap();
        let mut s = String::new();
        let _ = f.read_to_string(&mut s);
        let result: RunResult = json::decode(&s).unwrap();
        let length1 = result.strand1.length as i64;
        let length2 = result.strand2.length as i64;

        FlatPlotter {
            result: result,
            max_length: cmp::max(length1, length2) as f64,
            width: 1500.0,
            height: 900.0,
            svg: String::new(),
        }
    }

    fn plot_flat(&mut self, dx: f64, dy: f64) {
        self.svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                stroke='#222' stroke-width='10'/>
                "#,
                0,
                5,
                self.result.strand1.length as f64/self.max_length*self.width,
                5
                );
        self.svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                stroke='#222' stroke-width='10'/>
                "#,
                0,
                self.height-5.0,
                self.result.strand2.length as f64/self.max_length*self.width,
                self.height-5.0
                );

        for sd in &self.result.sds {
            let (mut left, mut right) = if !sd.reversed {
                (sd.left as f64, sd.right as f64)
            } else {
                let (a, b) = (sd.left as i64,
                              self.result.strand2.length as i64 - sd.right as i64 - sd.length as i64);
                (std::cmp::min(a, b) as f64, std::cmp::max(a, b) as f64)
            };

            if !(left >= 3031042417.0 || right >= 3031042417.0) || sd.reversed {
                continue
            }

            left = left/self.max_length * self.width;
            right = right/self.max_length * self.width;

            let color = Rgb::from_hsv(Hsv {
                hue: RgbHue::from_radians((left/self.width) as f64 * 2.0 * PI),
                saturation: 0.9,
                value: 0.9,
            });
            let color = format!("#{:x}{:x}{:x}",
                                (color.red * 255.0) as i8,
                                (color.green * 255.0) as i8,
                                (color.blue * 255.0) as i8);

            self.svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                fill='{}' fill-opacity='0.3' stroke='{}' stroke-opacity='0.9'/>
                "#,
                left,
                10,
                right,
                self.height-10.0,
                color,
                color);
        }
        // self.ticks(10000000, 40.0, 10, "#000000");
        // self.ticks(5000000, 20.0, 8, "#666666");
        // self.ticks(1000000, 10.0, 0, "#444444");
    }

    fn plot(&mut self) -> String {
        format!("<?xml version='1.0' encoding='iso-8859-1' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>{}</svg>",
                self.width,
                self.height,
                self.svg)
    }

    fn title(&mut self, x: i32, y: i32, text: &str, font_size: i32) {
        self.svg += &format!("<text x='{}' y='{}' font-family='Helvetica' \
                              font-size='{}'>{}</text>",
                             x,
                             y + font_size,
                             font_size,
                             text);
    }
}

fn plot_flat(filename: &str, out_filename: &str) {
    let width = 950.0;
    let height = 950.0;

    let mut plotter = FlatPlotter::new(filename, width);
    plotter.title(10, 10, filename, 20);
    plotter.plot_flat(0.0, 100.0);

    let mut f = File::create(out_filename).unwrap();
    f.write_all(plotter.plot().as_bytes()).unwrap();
}

fn main() {
    for argument in env::args().skip(1) {
        let out = &format!("{}.svg", &argument);
        plot_flat(&argument, out);
    }
}