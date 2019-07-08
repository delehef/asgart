extern crate rand;

use std::io::prelude::*;
use std::fs::File;
use std::f64::consts::PI;
use ::plot::*;
use separator::Separatable;

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
        let length = result.strand.length as f64;
        ChordPlotter {
            result,
            settings,
            length,
        }
    }

    fn plot(&self) -> Result<()> {
        let out_filename = format!("{}.svg", &self.settings.out_file);
        File::create(&out_filename)
            .and_then(|mut f| f.write_all(self.plot_chord().as_bytes()))
            .and_then(|_| { println!("Chord plot written to `{}`", &out_filename); Ok(()) })
            .chain_err(|| format!("Unable to write in `{}`", &out_filename))?;

        Ok(())
    }
}

impl ChordPlotter {
    fn angle(&self, x: f64) -> f64 {
        -x / self.length * 2.0 * PI
    }

    fn cartesian(&self, t: f64, r: f64) -> (f64, f64) {
        (CX + t.cos() * r, CY - t.sin() * r)
    }

    fn arc(&self, radius: f64, t1: f64, t2: f64) -> String {
        let (start_x, start_y) = self.cartesian(t1, radius);
        let (end_x, end_y) = self.cartesian(t2, radius);

        let large_flag = if t2 - t1 > PI/2.0 { 1 } else { 0 };
        let sweep_flag = if t2 - t1 > 0.0 { 0 } else { 1 };

        format!("M {} {} A {} {} {} {} {} {} {}",
                start_x, start_y,
                radius, radius, 0, large_flag, sweep_flag, end_x, end_y)
    }


    fn plot_chord(&self) -> String {
        let mut svg = String::new();
        svg += &format!("\n<g transform='translate({}, {})' >\n", 0, 0);

        for chr in &self.result.strand.map {
            let t1 = self.angle(chr.position as f64) - INTER_RING_SPACING;
            let t2 = self.angle(chr.position as f64 + chr.length as f64) + INTER_RING_SPACING;
            let tt = t1 + (t2-t1)/2.0;


            // Main chromosome bar
            svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='5' />\n",
                            self.arc(R + RING_WIDTH, t1, t2));
            if self.result.strand.map.len() > 1 {
                // Secondary chromosome bar
                svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='1.5' />\n",
                                self.arc(R + RING_WIDTH + OUT_CEILING * 0.7, t1, t2));
            }

            // Chromosomes labels
            let r = R + RING_WIDTH + RING_MARGIN;
            let (x, y) = self.cartesian(tt, r + if self.result.strand.map.len() > 1 { 65. } else { 20. });
            svg += &format!("<text x='{}' y='{}' font-family='Helvetica' font-size='8' fill='#333' transform='rotate({}, {}, {})'>\n{}\n</text>\n",
                            x, y,
                            -tt/(2.0*PI)*360.0 + 90.0, x, y,
                            str::replace(&chr.name, "chr", ""));
        }

        for family in &self.result.families {
            for sd in family {
                let (left, right) = (sd.global_left_position as i64, sd.global_right_position as i64);

                let t11 = self.angle(left as f64);
                let t12 = self.angle(left as f64 + sd.left_length as f64);
                let t1 = t11 + (t12 - t11)/2.0;

                let t21 = self.angle(right as f64);
                let t22 = self.angle(right as f64 + sd.right_length as f64);
                let t2 = t21 + (t22 - t21)/2.0;

                let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
                if width <= self.settings.min_thickness {width = self.settings.min_thickness};

                let color = if true {
                    // Direction-based color
                    if sd.reversed { &self.settings.color2 } else { &self.settings.color1 }
                } else {
                    // Position-based color
                    // let color = Rgb::from_hsv(Hsv {hue: RgbHue::from_radians(t1 + 1.5 * PI), saturation: 0.9, value: 0.9});
                    // format!("#{:x}{:x}{:x}", (color.red * 255.0) as u8, (color.green * 255.0) as u8, (color.blue * 255.0) as u8);
                    "#000000"
                };


                let ((x1, y1), (x2, y2), (cx, cy)) = if sd.chr_left != sd.chr_right || self.result.strand.map.len() == 1 {
                    let (x1, y1) = self.cartesian(t1, R);
                    let (x2, y2) = self.cartesian(t2, R);
                    let (cx, cy) = (CX, CY);
                    ((x1, y1), (x2, y2), (cx, cy))
                } else {
                    let tt = t1 + (t2-t1)/2.0;
                    let rin = R + RING_WIDTH + RING_MARGIN;
                    let rout = rin + OUT_CEILING;
                    let (x1, y1) = self.cartesian(t1, rin);
                    let (cx, cy) = self.cartesian(tt, rout);
                    let (x2, y2) = self.cartesian(t2, rin);
                    ((x1, y1), (x2, y2), (cx, cy))
                };


                let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
                svg += &format!(r#"
                                <path
                                d='{}' fill='none' stroke='{}' stroke-opacity='0.3' stroke-width='{}' class='sd'>
                                <title>{}</title>
                                </path>
                                "#,
                                path, color, width,
                                &format!("{}bp/{}bp\n{} → {}\n{} → {}",
                                         sd.left_length.separated_string(),
                                         sd.right_length.separated_string(),
                                         sd.global_left_position.separated_string(),
                                         (sd.global_left_position + sd.left_length).separated_string(),
                                         sd.global_right_position.separated_string(),
                                         (sd.global_right_position + sd.right_length).separated_string()
                                )
                );
            }
        }


        for features_family in &self.settings.feature_tracks {
            let color = format!("#{:2X}{:2X}{:2X}", rand::random::<i8>(), rand::random::<i8>(), rand::random::<i8>());
            for feature in features_family.iter() {
                for position in &feature.positions {
                    let (start, end) = match *position {
                        FeaturePosition::Relative { ref chr, start, length} => {
                            let chr = self.result.strand.find_chr(&chr).unwrap_or_else(|| panic!("Unable to find fragment `{}`", chr));
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
                    let font_size = 4.0;
                    svg += &format!("<polygon points='{},{} {},{} {},{} {},{}' style='fill:{};'/>\n",
                                    x0, y0,
                                    x1, y1,
                                    x2, y2,
                                    x3, y3,
                                    color
                    );
                    svg += &format!("<text x='{}' y='{}' font-family='Helvetica' font-size='{}' transform='rotate({}, {}, {})'>{}</text>",
                                    x3 + font_size, y3 + font_size, font_size,
                                    -t0/(2.0*PI)*360.0, x3, y3,
                                    feature.name);
                }
            }
        }

        svg += "</g>";
        format!("<?xml version='1.0' encoding='UTF-8'  standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'> \

                 <style type='text/css'> \
                 .sd:hover {{ stroke-opacity: 1.0; stroke: crimson; stroke-width: {}; }} \
                 </style> \

                 {} \
                 </svg>",
                    TOTAL_WIDTH, TOTAL_WIDTH,
                    2.0*self.settings.min_thickness,
                    svg
            )
    }
}
