use std::io::prelude::*;
use std::fs::File;
use std::f64::consts::PI;
use ::structs::*;
use ::plot::{Plotter, Settings};

const R: f64 = 200.0;

const RING_WIDTH: f64 = 5.0;
const RING_MARGIN: f64 = 5.0;

const OUT_CEILING: f64 = 20.0;
const INTER_RING_SPACING: f64 = 0.002;
const TOTAL_WIDTH: f64 = 2.0*(R + RING_MARGIN + RING_WIDTH + OUT_CEILING);

const CX: f64 = TOTAL_WIDTH/2.0;
const CY: f64 = TOTAL_WIDTH/2.0;



pub struct EyePlotter {
    settings: Settings,
    length: f64,
}

impl Plotter for EyePlotter {
    fn new(settings: Settings) -> EyePlotter {
        let length = settings.result.strand1.length as f64;
        EyePlotter {
            settings: settings,
            length: length,
        }
    }

    fn plot(self) {
        let mut f = File::create(&self.settings.out_file).expect(&format!("Unable to create `{}`", &self.settings.out_file));
        f.write_all(self.plot_chord().as_bytes()).unwrap();
    }
}

impl EyePlotter {
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

        // for chr in &self.settings.result.strand1.map {
        //     let t1 = self.angle(chr.position as f64) - INTER_RING_SPACING;
        //     let t2 = self.angle(chr.position as f64 + chr.length as f64) + INTER_RING_SPACING;
        //     let tt = t1 + (t2-t1)/2.0;

        //     // let tc1 = self.angle(chr.position as f64 + self.centromeres[&chr.name].0);
        //     // let tc2 = self.angle(chr.position as f64 + self.centromeres[&chr.name].1);

        //     let r = R + RING_WIDTH + RING_MARGIN;

        //     // Chromosomes bars
        //     // let (xbl, ybl) = self.cartesian(t1, R + RING_WIDTH);
        //     // let (xtl, ytl) = self.cartesian(t1, R + RING_WIDTH + 70.0);
        //     // let (xbr, ybr) = self.cartesian(t2, R + RING_WIDTH);
        //     // let (xtr, ytr) = self.cartesian(t2, R + RING_WIDTH + 70.0);
        //     // self.svg += &format!("<polygon points='{},{} {},{} {},{} {},{}' fill='#ccc' stroke-width='0'/>",
        //     //                      xbl, ybl,
        //     //                      xtl, ytl,
        //     //                      xtr, ytr,
        //     //                      xbr, ybr
        //     // );

        //     svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='5' />\n",
        //                          self.arc(R + RING_WIDTH, t1, t2)
        //     );
        //     // svg += &format!("<path d='{}' stroke='#afafaf' fill='none' stroke-width='5' />\n",
        //     //                      self.arc(R + RING_WIDTH, tc1, tc2)
        //     // );
        //     svg += &format!("<path d='{}' stroke='#ccc' fill='none' stroke-width='1.5' />\n",
        //                          self.arc(R + RING_WIDTH +58.0, t1, t2)
        //     );


        //     // Chromosomes labels
        //     let (x, y) = self.cartesian(tt, r + 65.0);
        //     svg += &format!("<text x='{}' y='{}' font-family='\"Helvetica\"' font-size='8' fill='#333' transform='rotate({}, {}, {})'>\n{}\n</text>\n",
        //                          x, y,
        //                          -tt/(2.0*PI)*360.0 + 90.0, x, y,
        //                          str::replace(&chr.name, "chr", ""));
        // }

        svg += &format!("<circle cx='{}' cy='{}' r='{}' fill='#001a44' filter='url(#blurMe)'/>", CX, CY, TOTAL_WIDTH/2.0 - 5.0);
        for sd in self.settings.result.sds.iter().filter(|&sd| self.inter_sd(sd)) {
            if sd.length < 5000 {
                continue;
            }
            let (left, right) = (sd.left as i64, sd.right as i64);

            let t11 = self.angle(left as f64);
            let t12 = self.angle(left as f64 + sd.length as f64);
            let mut t1 = t11 + (t12 - t11)/2.0;

            let t21 = self.angle(right as f64);
            let t22 = self.angle(right as f64 + sd.length as f64);
            let mut t2 = t21 + (t22 - t21)/2.0;

            let mut width = R * (2.0*(1.0 - (t12-t11).cos())).sqrt(); // Cf. Al-Kashi
            if width <= self.settings.thickness {width = self.settings.thickness};

            // let color = Rgb::from_hsv(Hsv {
            //     hue: RgbHue::from_radians(t1 + 1.5 * PI),
            //     saturation: 0.9,
            //     value: 0.9,
            // });

            let color = if sd.reversed {"#14889C"} else {"#ccc"};
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

            let f = border * (if ttt < PI/1.0{
                1.0
            } else {
                -1.0
            });
            let tt = t1 + (t2 - t1)/2.0;
            let cx = CX + tt.cos()*f;
            let cy = CY - tt.sin()*f;
            println!("t1\t {} t2\t {} - ttt\t {} sin(ttt)\t {} --> {}", t1, t2, ttt, ttt.sin(), f);
            println!("\t{} \t{} - \t{} \t{} --> \n",
                     t1/(2.0*PI)*360.0, t2/(2.0*PI)*360.0, ttt/(2.0*PI)*360.0, ttt.sin(),
            );

            // svg += &format!("<circle cx='{}' cy='{}' r='5' fill='red' />", cx, cy);
            let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
            svg += &format!("<path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>\n",
                            path, color, width);
        }
        svg += &format!("<circle cx='{}' cy='{}' r='{}' fill='#001a44' filter='url(#blurMe)'/>", CX, CY, 0.9*R);
        svg += &format!("<circle cx='{}' cy='{}' r='{}' fill='#222' filter='url(#blurMe)'/>", CX, CY, 0.45*R);
        for sd in self.settings.result.sds.iter().filter(|&sd| !self.inter_sd(sd)) {
            if sd.length < 3000 {
                continue;
            }
            let color = if sd.reversed {"#14889C"} else {"#ccc"};
            // let rin = R + RING_WIDTH + RING_MARGIN;
            // let rout = rin + OUT_CEILING;
            // let (x1, y1) = self.cartesian(t1, rin);
            // let (x2, y2) = self.cartesian(t2, rin);
            // let (cx, cy) = self.cartesian(tt, rout);
            // format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2)

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

            let (x1, y1) = self.cartesian(t1, 0.9*R);
            let (x2, y2) = self.cartesian(t2, 0.9*R);
            let (cx, cy) = (CX, CY); //self.cartesian(tt, rout);
            let path = format!("M {},{} Q {},{} {} {}", x1, y1, cx, cy, x2, y2);
            svg += &format!("<path d='{}' fill='none' stroke='{}' stroke-opacity='0.5' stroke-width='{}'/>\n",
                            path, color, width);
        }

        svg += &format!("<circle cx='{}' cy='{}' r='{}' fill='url(#rgrad)' filter='url(#blurMe)'/>", CX, CY, TOTAL_WIDTH/2.0);
        svg += "</g>";
        format!("<?xml version='1.0' encoding='iso-8859-1' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>\n{}\n \
                 <filter id='blurMe'> \
                 <feGaussianBlur in='SourceGraphic' stdDeviation='3' /> \
                 </filter> \
                 <defs> \
                    <radialGradient id='grad1' cx='50%' cy='50%' r='50%' fx='50%' fy='50%'> \
                        <stop offset='0%' style='stop-color:rgb(255,255,255); stop-opacity=0.1' /> \
                        <stop offset='100%' style='stop-color:rgb(0,0,255); stop-opacity=0.1' /> \
                    </radialGradient> \


    <radialGradient id='rgrad' cx='50%' cy='50%' r='75%' > \
    <stop offset='0%' style='stop-color:#001a44;stop-opacity:0' />\
<stop offset='38%' style='stop-color:#001a44;stop-opacity:0' />\
<stop offset='60%' style='stop-color:#001a44;stop-opacity:0.7' /> \
<stop offset='100%' style='stop-color:#001a44;stop-opacity:0.7' /> \
</radialGradient> \


    <mask id='fade' maskContentUnits='objectBoundingBox'> \
      <rect width='1' height='1' fill='url(#fadeGrad)'/> \
    </mask> \

                 </defs> \
                 </svg>",
                1.2*TOTAL_WIDTH,
                1.2*TOTAL_WIDTH,
                svg
        )
    }
}
