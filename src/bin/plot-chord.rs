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

const R: f64 = 400.0;
const CHORD_MAX_WIDTH :f64= 10.0;

const RING_WIDTH: f64 = 5.0;
const RING_MARGIN: f64 = 10.0;

const OUT_CEILING: f64 = 50.0;

const TOTAL_WIDTH: f64 = 2.0*(R + 2.0*RING_MARGIN + RING_WIDTH + OUT_CEILING);

const CX: f64 = TOTAL_WIDTH/2.0;
const CY: f64 = TOTAL_WIDTH/2.0;

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

    fn arc(&self, radius: f64, t1: f64, t2: f64) -> String {
        let (start_x, start_y) = self.cartesian(t1, radius);
        let (end_x, end_y) = self.cartesian(t2, radius);

        let flag = if t2 - t1 <= 180.0 {"0"} else {"1"};

        format!("M {} {} A {} {} {} {} {} {} {}", start_x, start_y, radius, radius, 0, flag, 1, end_x, end_y)
    }

    fn chr_left(&self, sd: &SD) -> isize {
        for (i, chr) in self.result.strand1.map.iter().enumerate() {
            if sd.left >= chr.position && sd.left <= chr.position + chr.length {
                return i as isize
            }
        }
        return -1;
    }

    fn chr_right(&self, sd: &SD) -> isize {
        for (i, chr) in self.result.strand2.map.iter().enumerate() {
            if sd.right >= chr.position && sd.right <= chr.position + chr.length {
                return i as isize
            }
        }
        return -1;
    }

    fn inter_sd(&self, sd: &SD) -> bool {
        self.chr_left(sd) != self.chr_right(sd)
    }

    fn plot_chord(&mut self, dx: f64, dy: f64) {
        self.svg += &format!("\n<g transform='translate({}, {})'>", dx, dy);
        {
            let radius = self.radius;
            self.ellipse(radius, radius, radius, radius);
        }

        for chr in &self.result.strand1.map {
            let t1 = self.angle(chr.position as f64) - 0.01;
            let t2 = self.angle(chr.position as f64+ chr.length as f64) + 0.01;

            self.svg += &format!("<path d='{}' stroke='#333' fill='none' stroke-width='5' />",
                                 self.arc(self.radius, t1, t2)
                                );

            let tt = t1 + (t2-t1)/2.0;
            let r = R + RING_WIDTH + RING_MARGIN;
            let (x, y) = self.cartesian(tt, self.radius + 30.0);
            println!("{}", tt/(2.0*PI)*360.0);
            self.svg += &format!("<text x='{}' y='{}' font-size='15' fill='#333' transform='rotate({}, {}, {})'>{}</text>", x, y, -tt/(2.0*PI)*360.0, x, y, chr.name);
        }

        for sd in &self.result.sds {
            let (left, right) = (sd.left as i64, sd.right as i64);

            let t1 = self.angle(left as f64);
            let t2 = self.angle(right as f64);
            let tt = t1 + (t2-t1)/2.0;

            let color = Rgb::from_hsv(Hsv {
                hue: RgbHue::from_radians(t1 + 1.5 * PI),
                saturation: 0.9,
                value: 0.9,
            });
            let color = format!("#{:x}{:x}{:x}",
                                (color.red * 255.0) as i8,
                                (color.green * 255.0) as i8,
                                (color.blue * 255.0) as i8);

            let path = if self.inter_sd(sd) {
                let x1 = CX + R * t1.cos();
                let y1 = CX + R * t1.sin();
                let x2 = CX + R * t2.cos();
                let y2 = CX + R * t2.sin();

                format!("M {},{} Q {},{} {} {}", x1, y1, CX, CY, x2, y2)
            } else {
                "".to_owned()
            };
            self.svg += &format!(" <path d='{}' fill='none' stroke='{}' stroke-opacity='0.9' stroke-width='1'/>", path, color);
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
