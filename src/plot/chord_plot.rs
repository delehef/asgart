extern crate rand;

use std::io::prelude::*;
use std::fs::File;
use std::collections::HashMap;
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
    centromeres: HashMap<String, (f64, f64)>,
}

impl Plotter for ChordPlotter {
    fn new(settings: Settings, result: RunResult) -> ChordPlotter {
        let mut _centromeres = HashMap::new();
        _centromeres.insert("chr1 ".to_owned(), (122_026_460.0, 125184587.0));
        _centromeres.insert("chr2 ".to_owned(), (92_188_146.0, 94090557.0));
        _centromeres.insert("chr3 ".to_owned(), (90_772_459.0, 93655574.0));
        _centromeres.insert("chr4 ".to_owned(), (49_708_101.0, 51743951.0));
        _centromeres.insert("chr5 ".to_owned(), (46_485_901.0, 50059807.0));
        _centromeres.insert("chr6 ".to_owned(), (58_553_889.0, 59829934.0));
        _centromeres.insert("chr7 ".to_owned(), (58_169_654.0, 60828234.0));
        _centromeres.insert("chr8 ".to_owned(), (44_033_745.0, 45877265.0));
        _centromeres.insert("chr9 ".to_owned(), (43_236_168.0, 45518558.0));
        _centromeres.insert("chr10 ".to_owned(), (39_686_383.0, 41593521.0));
        _centromeres.insert("chr11 ".to_owned(), (51_078_349.0, 54425074.0));
        _centromeres.insert("chr12 ".to_owned(), (34_769_408.0, 37185252.0));
        _centromeres.insert("chr13 ".to_owned(), (16_000_001.0, 18051248.0));
        _centromeres.insert("chr14 ".to_owned(), (16_000_001.0, 18173523.0));
        _centromeres.insert("chr15 ".to_owned(), (17_000_001.0, 19725254.0));
        _centromeres.insert("chr16 ".to_owned(), (36_311_159.0, 38280682.0));
        _centromeres.insert("chr17 ".to_owned(), (22_813_680.0, 26885980.0));
        _centromeres.insert("chr18 ".to_owned(), (15_460_900.0, 20861206.0));
        _centromeres.insert("chr19 ".to_owned(), (24_498_981.0, 27190874.0));
        _centromeres.insert("chr20 ".to_owned(), (26_436_233.0, 30038348.0));
        _centromeres.insert("chr21 ".to_owned(), (10_864_561.0, 12915808.0));
        _centromeres.insert("chr22 ".to_owned(), (12_954_789.0, 15054318.0));
        _centromeres.insert("chrX ".to_owned(), (58_605_580.0, 62412542.0));
        _centromeres.insert("chrY ".to_owned(), (10_316_745.0, 10544039.0));


        let length = result.strand1.length as f64;
        ChordPlotter {
            result: result,
            settings: settings,

            length: length,
            centromeres: _centromeres,
        }
    }

    fn plot(self) {
        let mut f = File::create(&self.settings.out_file).expect(&format!("Unable to create `{}`", &self.settings.out_file));
        f.write_all(self.plot_chord().as_bytes()).unwrap();
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

    fn intra_sd(&self, sd: &SD) -> bool {
        !self.inter_sd(sd)
    }

    fn plot_chord(self) -> String {
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


            if false { // Plot centromeres ?
                // Start/end angles of centromeres
                let tc1 = self.angle(chr.position as f64 + self.centromeres[&chr.name].0);
                let tc2 = self.angle(chr.position as f64 + self.centromeres[&chr.name].1);
                svg += &format!("<path d='{}' stroke='#afafaf' fill='none' stroke-width='5' />\n",
                                self.arc(R + RING_WIDTH, tc1, tc2));
            }


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
                let (left, right) = (sd.left as i64, sd.right as i64);

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
                let (left, right) = (sd.left as i64, sd.right as i64);

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
