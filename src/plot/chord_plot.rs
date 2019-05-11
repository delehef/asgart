extern crate rand;

use std::io::prelude::*;
use std::fs::File;
use std::f64::consts::PI;
use ::plot::*;

const R: f64 = 200.0;

const RING_WIDTH: f64 = 5.0;
const RING_MARGIN: f64 = 10.0;

const OUT_CEILING: f64 = R/2.0;
const INTER_RING_SPACING: f64 = 0.002;
const TOTAL_WIDTH: f64 = 2.5*(R + RING_MARGIN + RING_WIDTH + OUT_CEILING);

const CX: f64 = TOTAL_WIDTH/2.0;
const CY: f64 = TOTAL_WIDTH/2.0;



pub struct ChordPlotter {
    result: RunResult,
    settings: Settings,

    length: f64,
}

impl Plotter for ChordPlotter {
    fn new(settings: Settings, result: RunResult) -> ChordPlotter {
        let length = result.strand1.length as f64;
        ChordPlotter {
            result: result,
            settings: settings,

            length: length,
        }
    }

    fn plot(&self) -> Result<()> {
        let out_filename = format!("{}.svg", &self.settings.out_file);
        File::create(&out_filename)
            .and_then(|mut f| f.write_all(self.plot_chord().as_bytes()).into())
            .and_then(|_| Ok(println!("Chord plot written to `{}`", &out_filename)))
            .chain_err(|| format!("Unable to write in `{}`", &out_filename))?;

        Ok(())
    }
}

impl ChordPlotter {
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

    fn inter_sd(&self, sd: &SD) -> bool {sd.chr_left != sd.chr_right}
    fn intra_sd(&self, sd: &SD) -> bool {!self.inter_sd(sd)}

    fn plot_chord(&self) -> String {
        let mut svg = String::new();
        svg += &format!("\n<g transform='translate({}, {})' >\n", 0, 0);

        for chr in &self.result.strand1.map {
            let t1 = self.angle(chr.position as f64) - INTER_RING_SPACING;
            let t2 = self.angle(chr.position as f64 + chr.length as f64) + INTER_RING_SPACING;
            let tt = t1 + (t2-t1)/2.0;


            // Main chromosome bar
            svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='5' />\n",
                            self.arc(R + RING_WIDTH, t1, t2));
            // Secondary chromosome bar
            svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='1.5' />\n",
                            self.arc(R + RING_WIDTH + OUT_CEILING * 0.7, t1, t2));

            // Chromosomes labels
            let r = R + RING_WIDTH + RING_MARGIN;
            let (x, y) = self.cartesian(tt, r + 65.0);
            svg += &format!("<text x='{}' y='{}' font-family='\"Helvetica\"' font-size='8' fill='#333' transform='rotate({}, {}, {})'>\n{}\n</text>\n",
                            x, y,
                            -tt/(2.0*PI)*360.0 + 90.0, x, y,
                            str::replace(&chr.name, "chr", ""));
        }

        for sd in self.result.sds
            .iter()
            .filter(|&sd| self.inter_sd(sd)) {
                let (left, right) = (sd.global_left_position as i64, sd.global_right_position as i64);

                let t11 = self.angle(left as f64);
                let t12 = self.angle(left as f64 + sd.length as f64);
                let mut t1 = t11 + (t12 - t11)/2.0;

                let t21 = self.angle(right as f64);
                let t22 = self.angle(right as f64 + sd.length as f64);
                let mut t2 = t21 + (t22 - t21)/2.0;

                let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
                if width <= self.settings.thickness {width = self.settings.thickness};

                let color = if true {
                    // Direction-based color
                    if sd.reversed { &self.settings.color2 } else { &self.settings.color1 }
                } else {
                    // Position-based color
                    // let color = Rgb::from_hsv(Hsv {
                    //     hue: RgbHue::from_radians(t1 + 1.5 * PI),
                    //     saturation: 0.9,
                    //     value: 0.9,
                    // });
                    // format!("#{:x}{:x}{:x}",
                    //         (color.red * 255.0) as u8,
                    //         (color.green * 255.0) as u8,
                    //         (color.blue * 255.0) as u8);
                    "#000000"
                };



                let (x1, y1) = self.cartesian(t1, R);
                let (x2, y2) = self.cartesian(t2, R);

                while t2 < 0.0 {t2 += 2.0*PI;}
                while t2 > 2.0*PI {t2 -= 2.0*PI;}
                while t1 > 2.0*PI {t1 -= 2.0*PI;}

                let cx = CX;
                let cy = CY;

                let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
                svg += &format!("<path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>\n",
                                path, color, width);
            }

        for sd in self.result.sds
            .iter()
            .filter(|&sd| self.intra_sd(sd)) {
                let (left, right) = (sd.global_left_position as i64, sd.global_left_position as i64);

                let t11 = self.angle(left as f64);
                let t12 = self.angle(left as f64 + sd.length as f64);
                let t1 = t11 + (t12 - t11)/2.0;

                let t21 = self.angle(right as f64);
                let t22 = self.angle(right as f64 + sd.length as f64);
                let t2 = t21 + (t22 - t21)/2.0;

                let tt = t1 + (t2-t1)/2.0;
                let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
                if width <= self.settings.thickness {width = self.settings.thickness};

                let color = if true {
                    // Direction-based color
                    if sd.reversed { &self.settings.color2 } else { &self.settings.color1 }
                } else {
                    // Position-based color
                    // let color = Rgb::from_hsv(Hsv {
                    //     hue: RgbHue::from_radians(t1 + 1.5 * PI),
                    //     saturation: 0.9,
                    //     value: 0.9,
                    // });
                    // format!("#{:x}{:x}{:x}",
                    //         (color.red * 255.0) as u8,
                    //         (color.green * 255.0) as u8,
                    //         (color.blue * 255.0) as u8);
                    "#000000"
                };
                let rin = R + RING_WIDTH + RING_MARGIN;
                let rout = rin + OUT_CEILING;
                let (x1, y1) = self.cartesian(t1, rin);
                let (x2, y2) = self.cartesian(t2, rin);
                let (cx, cy) = self.cartesian(tt, rout);
                let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
                svg += &format!("<path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>\n",
                                path, color, width);
            }

        for features_family in &self.settings.feature_tracks {
            let color = format!("#{:2X}{:2X}{:2X}", rand::random::<i8>(), rand::random::<i8>(), rand::random::<i8>());
            for feature in features_family.iter() {
                for position in &feature.positions {
                    let (start, end) = match *position {
                        FeaturePosition::Relative { ref chr, start, length} => {
                            let chr = self.result.strand1.find_chr(&chr);
                            (chr.position + start, chr.position + start + length)
                        }
                        FeaturePosition::Absolute { start, length }         => { (start, start + length) }
                    };
                    let t1 = self.angle(start as f64);
                    let t2 = self.angle(end as f64);
                    let t0 = t1 + (t2 - t1)/2.0;

                    let (x0, y0) = self.cartesian(t0 - 0.02, R - 5.0);
                    let (x1, y1) = self.cartesian(t1, R);
                    let (x2, y2) = self.cartesian(t2, R);
                    let (x3, y3) = self.cartesian(t0 + 0.02, R - 5.0);
                    let font_size = 6.0;
                    svg += &format!("<polygon points='{},{} {},{} {},{} {},{}' style='fill:{};'/>\n",
                                    x0, y0,
                                    x1, y1,
                                    x2, y2,
                                    x3, y3,
                                    color
                    );
                    svg += &format!("<text x='{}' y='{}' font-family='Helvetica' font-size='{}'>{}</text>",
                                    x3, y3 + font_size,
                                    font_size, feature.name);
                }
            }
        }

        svg += "</g>";
        format!("<?xml version='1.0' encoding='iso-8859-1' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>\n{}\n \
                 </svg>",
                TOTAL_WIDTH,
                TOTAL_WIDTH,
                svg
        )
    }
}
