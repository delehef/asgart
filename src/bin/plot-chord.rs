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
const CHORD_MIN_WIDTH: f64 = 0.2;

const RING_WIDTH: f64 = 5.0;
const RING_MARGIN: f64 = 10.0;

const OUT_CEILING: f64 = 100.0;

const TOTAL_WIDTH: f64 = 3.0*(R + 2.0*RING_MARGIN + RING_WIDTH + OUT_CEILING);

const CX: f64 = TOTAL_WIDTH/2.0;
const CY: f64 = TOTAL_WIDTH/2.0;

use asgart::structs::*;



struct ChordPlotter {
    result: RunResult,
    length: f64,
    radius: f64,
    svg: String,
}

impl ChordPlotter {
    fn new(json_file: &str, radius: f64) -> ChordPlotter {
        let mut f = File::open(json_file).unwrap();
        let mut s = String::new();
        let _ = f.read_to_string(&mut s);
        let result: RunResult = json::decode(&s).unwrap();

        ChordPlotter {
            length: result.strand1.length as f64,
            result: result,

            radius: radius,
            svg: String::new(),
        }
    }

    fn angle(&self, x: f64) -> f64 {
         PI/2.0 - x / self.length * 2.0 * PI
    }

    fn cartesian(&self, t: f64, r: f64) -> (f64, f64) {
        (CX + t.cos() * r, CY - t.sin() * r)
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
        -1
    }

    fn chr_right(&self, sd: &SD) -> isize {
        for (i, chr) in self.result.strand2.map.iter().enumerate() {
            if sd.right >= chr.position && sd.right <= chr.position + chr.length {
                return i as isize
            }
        }
        -1
    }

    fn inter_sd(&self, sd: &SD) -> bool {
        self.chr_left(sd) != self.chr_right(sd)
    }

    fn plot_chord(&mut self, dx: f64, dy: f64) {
        self.svg += &format!("\n<g transform='translate({}, {})'>", dx, dy);

        for chr in &self.result.strand1.map {
            let t1 = self.angle(chr.position as f64) - 0.005;
            let t2 = self.angle(chr.position as f64+ chr.length as f64) + 0.005;

            self.svg += &format!("<path d='{}' stroke='#333' fill='none' stroke-width='5' />",
                                 self.arc(self.radius, t1, t2)
                                );

            let tt = t1 + (t2-t1)/2.0;
            let r = R + RING_WIDTH + RING_MARGIN;
            let (x, y) = self.cartesian(tt, r + 30.0);
            self.svg += &format!("<text x='{}' y='{}' font-size='15' fill='#333' transform='rotate({}, {}, {})'>{}</text>", x, y, -tt/(2.0*PI)*360.0, x, y, chr.name);
        }

        for sd in &self.result.sds {
            if (sd.identity > 0.0 && sd.identity < 99.0) || sd.length < 10000 || !sd.reversed {
                continue;
            }
            let (left, right) = (sd.left as i64, sd.right as i64);

            let t11 = self.angle(left as f64);
            let t12 = self.angle(left as f64 + sd.length as f64);
            let t1 = t11 + (t12 - t11)/2.0;

            let t21 = self.angle(right as f64);
            let t22 = self.angle(right as f64 + sd.length as f64);
            let t2 = t21 + (t22 - t21)/2.0;

            let tt = t1 + (t2-t1)/2.0;
            let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
            if width <= CHORD_MIN_WIDTH {width = CHORD_MIN_WIDTH};

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
                let (x1, y1) = self.cartesian(t1, R);
                let (x2, y2) = self.cartesian(t2, R);

                format!("M {},{} Q {},{} {} {}", x1, y1, CX, CY, x2, y2)
            } else {
                let rin = R + RING_WIDTH + RING_MARGIN;
                let rout = rin + OUT_CEILING;
                let (x1, y1) = self.cartesian(t1, rin);
                let (x2, y2) = self.cartesian(t2, rin);
                let (cx, cy) = self.cartesian(tt, rout);
                format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2)
            };
            self.svg += &format!(" <path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>", path, color, width);
        }
        // self.ticks();
        self.svg += "</g>";
    }

    fn plot(&mut self) -> String {
        format!("<?xml version='1.0' encoding='iso-8859-1' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>{}</svg>",
                 TOTAL_WIDTH,
                 TOTAL_WIDTH,
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

fn plot_chord(filename: &str, out_filename: &str) {
    let radius = 400.0;

    let mut plotter = ChordPlotter::new(filename, radius);
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
