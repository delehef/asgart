use separator::Separatable;

use std::cmp;
use std::io::prelude::*;
use std::fs::File;
use ::plot::{Plotter, Settings};

const CHR_WIDTH: f64 = 4.0;

pub struct FlatPlotter {
    settings: Settings,
    max_length: f64,
    width: f64,
    height: f64,
}

impl Plotter for FlatPlotter {
    fn new(settings: Settings) -> FlatPlotter {
        let length1 = settings.result.strand1.length as i64;
        let length2 = settings.result.strand2.length as i64;

        FlatPlotter {
            settings: settings,
            max_length: cmp::max(length1, length2) as f64,
            width: 1500.0,
            height: 230.0,
        }
    }

    fn plot(self) {
        let mut f = File::create(&self.settings.out_file).unwrap();
        f.write_all(self.plot_flat().as_bytes()).unwrap();
    }
}

impl FlatPlotter {
    fn plot_flat(self) -> String {
        let mut svg = String::new();
        //
        // Chromosomes
        //
        svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                stroke='#ccc' stroke-width='{}'/>
                "#,
                             0,
                             CHR_WIDTH/2.0,
                             self.settings.result.strand1.length as f64/self.max_length*self.width,
                             CHR_WIDTH/2.0,
                             CHR_WIDTH
        );
        svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                stroke='#ccc' stroke-width='{}'/>
                "#,
                             0,
                             self.height-CHR_WIDTH/2.0,
                             self.settings.result.strand2.length as f64/self.max_length*self.width,
                             self.height-CHR_WIDTH/2.0,
                             CHR_WIDTH
        );
        let centromere_start = 10316945.0;
        let centromere_end = 10544039.0;
        svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                stroke='#afafaf' stroke-width='{}'/>
                "#,
                             centromere_start as f64/self.max_length*self.width,
                             CHR_WIDTH/2.0,
                             centromere_end as f64/self.max_length*self.width,
                             CHR_WIDTH/2.0,
                             CHR_WIDTH
        );
        svg += &format!(r#"
                <line
                x1='{}' y1='{}' x2='{}' y2='{}'
                stroke='#afafaf' stroke-width='{}'/>
                "#,
                             centromere_start as f64/self.max_length*self.width,
                             self.height-CHR_WIDTH/2.0,
                             centromere_end as f64/self.max_length*self.width,
                             self.height-CHR_WIDTH/2.0,
                             CHR_WIDTH
        );



        //
        // Ticks
        //
        for i in 0..self.max_length as i64 {
            if i % 1000000 == 0 {
                let height = if i % 10000000 == 0 {
                    self.height + 15.0
                } else if i % 5000000 == 0 {
                    self.height + 10.0
                } else {
                    self.height + 5.0
                };
                let x = i as f64/self.max_length*self.width;

                svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#898989' stroke-width='1'/>",
                                    x, self.height,
                                    x, height
                );

                if i % 10000000 == 0 {
                    svg += &format!("<text x='{}' y='{}' font-family='Helvetica' font-size='8'>{}Mb</text>",
                                         x+4.0,
                                         height,
                                         i / 1000000);
                }
            }
        }

        for sd in self.settings.result.sds.iter().filter(|&sd| sd.length >= self.settings.min_length) {
            let left1 = (sd.left as f64)/self.max_length * self.width;
            let left2 = (sd.left as f64 + sd.length as f64)/self.max_length * self.width;
            let right1 = (sd.right as f64)/self.max_length * self.width;
            let right2 = (sd.right as f64+ sd.length as f64)/self.max_length * self.width;

            let color = if sd.reversed {"#00b2ae"} else {"#ff5b00"};

            svg += &format!(r#"
                            <polygon
                            points='{},{} {},{} {},{} {},{}'
                            fill='{}' fill-opacity='0.5' stroke='{}' stroke-opacity='0.9'
                            stroke-width='0'>
                            >
                            <title>{}</title>
                            </polygon>
                            "#,
                            left1, CHR_WIDTH,
                            left2, CHR_WIDTH,
                            right2, self.height - CHR_WIDTH,
                            right1, self.height - CHR_WIDTH,
                            color,
                            color,
                            &format!("{}bp\n{} → {}\n{} → {}",
                                     sd.length.separated_string(),
                                     sd.left.separated_string(),
                                     (sd.left + sd.length).separated_string(),
                                     sd.right.separated_string(),
                                     (sd.right + sd.length).separated_string()
                            )
            );
        }
        format!("<?xml version='1.0' encoding='UTF-8' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>{}</svg>",
                self.width,
                self.height+30.0,
                svg)
    }
}
