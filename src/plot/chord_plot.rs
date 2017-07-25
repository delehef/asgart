use std::io::prelude::*;
use std::fs::File;
use std::path::Path;
use std::collections::HashMap;
use std::f64::consts::PI;
use ::structs::*;
use ::plot::{Plotter, Settings};

const R: f64 = 200.0;

const RING_WIDTH: f64 = 5.0;
const RING_MARGIN: f64 = 10.0;

const OUT_CEILING: f64 = R/2.0;
const INTER_RING_SPACING: f64 = 0.002;
const TOTAL_WIDTH: f64 = 2.5*(R + RING_MARGIN + RING_WIDTH + OUT_CEILING);

const CX: f64 = TOTAL_WIDTH/2.0;
const CY: f64 = TOTAL_WIDTH/2.0;



pub struct ChordPlotter {
    settings: Settings,
    length: f64,
    centromeres: HashMap<String, (f64, f64)>,
}

impl Plotter for ChordPlotter {
    fn new(settings: Settings) -> ChordPlotter {
        let mut centromeres = HashMap::new();
        centromeres.insert("chr1 ".to_owned(), (122026460.0, 125184587.0));
        centromeres.insert("chr2 ".to_owned(), (92188146.0, 94090557.0));
        centromeres.insert("chr3 ".to_owned(), (90772459.0, 93655574.0));
        centromeres.insert("chr4 ".to_owned(), (49708101.0, 51743951.0));
        centromeres.insert("chr5 ".to_owned(), (46485901.0, 50059807.0));
        centromeres.insert("chr6 ".to_owned(), (58553889.0, 59829934.0));
        centromeres.insert("chr7 ".to_owned(), (58169654.0, 60828234.0));
        centromeres.insert("chr8 ".to_owned(), (44033745.0, 45877265.0));
        centromeres.insert("chr9 ".to_owned(), (43236168.0, 45518558.0));
        centromeres.insert("chr10 ".to_owned(), (39686383.0, 41593521.0));
        centromeres.insert("chr11 ".to_owned(), (51078349.0, 54425074.0));
        centromeres.insert("chr12 ".to_owned(), (34769408.0, 37185252.0));
        centromeres.insert("chr13 ".to_owned(), (16000001.0, 18051248.0));
        centromeres.insert("chr14 ".to_owned(), (16000001.0, 18173523.0));
        centromeres.insert("chr15 ".to_owned(), (17000001.0, 19725254.0));
        centromeres.insert("chr16 ".to_owned(), (36311159.0, 38280682.0));
        centromeres.insert("chr17 ".to_owned(), (22813680.0, 26885980.0));
        centromeres.insert("chr18 ".to_owned(), (15460900.0, 20861206.0));
        centromeres.insert("chr19 ".to_owned(), (24498981.0, 27190874.0));
        centromeres.insert("chr20 ".to_owned(), (26436233.0, 30038348.0));
        centromeres.insert("chr21 ".to_owned(), (10864561.0, 12915808.0));
        centromeres.insert("chr22 ".to_owned(), (12954789.0, 15054318.0));
        centromeres.insert("chrX ".to_owned(), (58605580.0, 62412542.0));
        centromeres.insert("chrY ".to_owned(), (10316745.0, 10544039.0));


        let length = settings.result.strand1.length as f64;
        ChordPlotter {
            settings: settings,
            length: length,

            centromeres: centromeres.clone(),
        }
    }

    fn plot(self) {
        let mut f = File::create(&self.settings.out_file).expect(&format!("Unable to create `{}`", &self.settings.out_file));
        f.write_all(self.plot_chord().as_bytes()).unwrap();
    }
}

impl ChordPlotter {
    fn title(self, x: i32, y: i32, text: &str, font_size: i32) -> String {
        format!("<text x='{}' y='{}' font-family='Helvetica' \
                          font-size='{}'>{}</text>",
                          x,
                          y + font_size,
                          font_size,
                          text)
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
        for (i, chr) in self.settings.result.strand1.map.iter().enumerate() {
            if sd.left >= chr.position && sd.left <= chr.position + chr.length {
                return i as isize
            }
        }
        -1
    }

    fn chr_right(&self, sd: &SD) -> isize {
        for (i, chr) in self.settings.result.strand2.map.iter().enumerate() {
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

        for chr in &self.settings.result.strand1.map {
            let t1 = self.angle(chr.position as f64) - INTER_RING_SPACING;
            let t2 = self.angle(chr.position as f64 + chr.length as f64) + INTER_RING_SPACING;
            let tt = t1 + (t2-t1)/2.0;

            // let tc1 = self.angle(chr.position as f64 + self.centromeres[&chr.name].0);
            // let tc2 = self.angle(chr.position as f64 + self.centromeres[&chr.name].1);

            let r = R + RING_WIDTH + RING_MARGIN;

            // Chromosomes bars
            // let (xbl, ybl) = self.cartesian(t1, R + RING_WIDTH);
            // let (xtl, ytl) = self.cartesian(t1, R + RING_WIDTH + 70.0);
            // let (xbr, ybr) = self.cartesian(t2, R + RING_WIDTH);
            // let (xtr, ytr) = self.cartesian(t2, R + RING_WIDTH + 70.0);
            // self.svg += &format!("<polygon points='{},{} {},{} {},{} {},{}' fill='#ccc' stroke-width='0'/>",
            //                      xbl, ybl,
            //                      xtl, ytl,
            //                      xtr, ytr,
            //                      xbr, ybr
            // );

            svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='5' />\n",
                                 self.arc(R + RING_WIDTH, t1, t2)
            );
            // svg += &format!("<path d='{}' stroke='#afafaf' fill='none' stroke-width='5' />\n",
            //                      self.arc(R + RING_WIDTH, tc1, tc2)
            // );
            svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='1.5' />\n",
                                 self.arc(R + RING_WIDTH + OUT_CEILING * 0.7, t1, t2)
            );


            // Chromosomes labels
            let (x, y) = self.cartesian(tt, r + 65.0);
            svg += &format!("<text x='{}' y='{}' font-family='\"Helvetica\"' font-size='8' fill='#333' transform='rotate({}, {}, {})'>\n{}\n</text>\n",
                                 x, y,
                                 -tt/(2.0*PI)*360.0 + 90.0, x, y,
                                 str::replace(&chr.name, "chr", ""));
        }

        for sd in self.settings.result.sds.iter().filter(|&sd| self.inter_sd(sd) && sd.length >= self.settings.min_length) {
            let (left, right) = (sd.left as i64, sd.right as i64);

            let t11 = self.angle(left as f64);
            let t12 = self.angle(left as f64 + sd.length as f64);
            let mut t1 = t11 + (t12 - t11)/2.0;

            let t21 = self.angle(right as f64);
            let t22 = self.angle(right as f64 + sd.length as f64);
            let mut t2 = t21 + (t22 - t21)/2.0;

            let tt = t1 + (t2-t1)/2.0;
            let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
            if width <= self.settings.thickness {width = self.settings.thickness};

            // let color = Rgb::from_hsv(Hsv {
            //     hue: RgbHue::from_radians(t1 + 1.5 * PI),
            //     saturation: 0.9,
            //     value: 0.9,
            // });

            let color = if sd.reversed {"#00b2ae"} else {"#ff5b00"};
            // let color = format!("#{:x}{:x}{:x}",
            //                     (color.red * 255.0) as u8,
            //                     (color.green * 255.0) as u8,
            //                     (color.blue * 255.0) as u8);

            let border = 1.4*R;
            let (x1, y1) = self.cartesian(t1, R);
            let (x2, y2) = self.cartesian(t2, R);
            
            while t2 < 0.0 {t2 += 2.0*PI;}
            while t2 > 2.0*PI {t2 -= 2.0*PI;}

            while t1 > 2.0*PI {t1 -= 2.0*PI;}

            let (t1, t2) = if t2 < t1 {(t2, t1)} else {(t1, t2)};

            let mut ttt = t2 - t1;
            
            let tt = t1 + (t2 - t1)/2.0;
            let cx = CX;
            let cy = CY;

            let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
            svg += &format!("<path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>\n",
                            path, color, width);
        }

        for sd in self.settings.result.sds.iter().filter(|&sd| self.intra_sd(sd) && sd.length >= self.settings.min_length) {
            let (left, right) = (sd.left as i64, sd.right as i64);

            let t11 = self.angle(left as f64);
            let t12 = self.angle(left as f64 + sd.length as f64);
            let mut t1 = t11 + (t12 - t11)/2.0;

            let t21 = self.angle(right as f64);
            let t22 = self.angle(right as f64 + sd.length as f64);
            let mut t2 = t21 + (t22 - t21)/2.0;

            let tt = t1 + (t2-t1)/2.0;
            let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
            if width <= self.settings.thickness {width = self.settings.thickness};

            // let color = Rgb::from_hsv(Hsv {
            //     hue: RgbHue::from_radians(t1 + 1.5 * PI),
            //     saturation: 0.9,
            //     value: 0.9,
            // });

            let color = if sd.reversed {"#00b2ae"} else {"#ff5b00"};
            let rin = R + RING_WIDTH + RING_MARGIN;
            let rout = rin + OUT_CEILING;
            let (x1, y1) = self.cartesian(t1, rin);
            let (x2, y2) = self.cartesian(t2, rin);
            let (cx, cy) = self.cartesian(tt, rout);
            let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
            svg += &format!("<path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>\n",
                            path, color, width);
        }

        svg += "</g>";
        let title = format!("{} (>{}Kbp)",
            Path::new(&self.settings.out_file).file_stem().unwrap().to_str().unwrap().to_owned(),
            self.settings.min_length/1000
            );

        svg += &(self.title(10, 10, &title, 20));
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
