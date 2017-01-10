#![feature(step_by)]

extern crate asgart;
extern crate rustc_serialize;
extern crate palette;

use std::env;
use std::io::prelude::*;
use std::fs::File;
use rustc_serialize::json;
use std::f64::consts::PI;
use palette::{Rgb, Hsv, RgbHue, FromColor};

use asgart::structs::*;

struct ChordPlotter {
    result: RunResult,
    length: f64,
    width: f64,
    height: f64,
    radius: f64,
    svg: String,
}

impl ChordPlotter {
    fn new(json_file: &str, width: f64, height: f64, radius: f64) -> ChordPlotter {
        let mut f = File::open(json_file).unwrap();
        let mut s = String::new();
        let _ = f.read_to_string(&mut s);
        let result: RunResult = json::decode(&s).unwrap();

        ChordPlotter {
            length: result.strand1.length as f64,
            result: result,
            width: width,
            height: height,
            radius: radius,
            svg: String::new(),
        }
    }

    fn ticks(&mut self, step: i64, size: f64, font_size: i32, color: &str) {
        for l in (0..self.length as i64).step_by(step) {
            let t = PI / 2.0 - (l as f64) / self.length * 2.0 * PI;
            let x1 = self.radius + ((self.radius + 10.0) * t.cos());
            let x2 = self.radius + ((self.radius + 10.0 + size) * t.cos());
            let y1 = self.radius - ((self.radius + 10.0) * t.sin());
            let y2 = self.radius - ((self.radius + 10.0 + size) * t.sin());

            if font_size > 0 {
                self.svg += &format!("<text x='{}' y='{}' font-family='Helvetica' \
                                      font-size='{}'>{}M</text>",
                                     x2,
                                     y2 as i32 - font_size,
                                     font_size,
                                     l / 1000000);
            }

            self.svg += &format!("<path d='M{} {} L {} {}' stroke='{}'/>",
                                 x1,
                                 y1,
                                 x2,
                                 y2,
                                 color);
        }
    }

    fn angle(&self, x: f64) -> f64 {
        PI / 2.0 - x / self.length * 2.0 * PI
    }

    fn cartesian(&self, t: f64, r: f64) -> (f64, f64) {
        (self.radius + t.cos() * r, self.radius - t.sin() * r)
    }

    fn plot_chord(&mut self, dx: f64, dy: f64) {
        self.svg += &format!("\n<g transform='translate({}, {})'>", dx, dy);
        {
            let radius = self.radius;
            self.ellipse(radius, radius, radius, radius);
        }

        for sd in &self.result.sds {
            let (left, right) = if !sd.reversed {
                (sd.left as i64, sd.right as i64)
            } else {
                let (a, b) = (sd.left as i64,
                              self.length as i64 - sd.right as i64 - sd.length as i64);
                (std::cmp::min(a, b), std::cmp::max(a, b))
            };

            let t_left1 = self.angle(left as f64);
            let t_right2 = self.angle(right as f64);
            let t_left2 = self.angle((left + sd.length as i64) as f64);
            let t_right1 = self.angle((right + sd.length as i64) as f64);

            let pleft1 = self.cartesian(t_left1, self.radius);
            let pleft2 = self.cartesian(t_left2, self.radius);
            let pright1 = self.cartesian(t_right1, self.radius);
            let pright2 = self.cartesian(t_right2, self.radius);

            // let gradient = Gradient::new(vec![Hsv::new(0.0, 1.0, 1.0), Hsv::new(1.0, 1.0, 1.0)])
            let color = Rgb::from_hsv(Hsv {
                hue: RgbHue::from_radians(t_left1 + 1.5 * PI),
                saturation: 0.9,
                value: 0.9,
            });
            let color = format!("#{:x}{:x}{:x}",
                                (color.red * 255.0) as i8,
                                (color.green * 255.0) as i8,
                                (color.blue * 255.0) as i8);

            self.svg += &format!(r#"
                <path d='
                M{},{}
                Q{},{}
                ,{} {}

                L{},{}

                Q{},{}
                ,{},{}

                L{},{}
                '
                fill='{}' fill-opacity='0.3'
                stroke='{}' stroke-opacity='0.9'/>
                "#,
                                 pleft1.0,
                                 pleft1.1,
                                 self.radius,
                                 self.radius,
                                 pright1.0,
                                 pright1.1,
                                 pright2.0,
                                 pright2.1,
                                 self.radius,
                                 self.radius,
                                 pleft2.0,
                                 pleft2.1,
                                 pleft1.0,
                                 pleft1.1,
                                 color,
                                 color);
        }
        self.ticks(10000000, 40.0, 10, "#000000");
        self.ticks(5000000, 20.0, 8, "#666666");
        self.ticks(1000000, 10.0, 0, "#444444");
        self.svg += "</g>";
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

    fn ellipse(&mut self, cx: f64, cy: f64, rx: f64, ry: f64) {
        self.svg += &format!("<ellipse cx='{}' cy='{}' rx='{}' ry='{}' style='stroke:black; \
                              fill: none'/>",
                             cx,
                             cy,
                             rx,
                             ry);
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

fn plot_chord(filename: &str, out_filename: &str) {
    let width = 950.0;
    let height = 950.0;
    let radius = 400.0;

    let mut plotter = ChordPlotter::new(filename, width, height, radius);
    plotter.title(10, 10, filename, 20);
    plotter.plot_chord(0.0, 100.0);

    let mut f = File::create(out_filename).unwrap();
    f.write_all(plotter.plot().as_bytes()).unwrap();
}

fn main() {
    for argument in env::args().skip(1) {
        let out = &format!("{}.svg", &argument);
        plot_chord(&argument, out);
    }
}
