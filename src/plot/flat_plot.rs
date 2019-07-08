extern crate rand;

use separator::Separatable;

use std::io::prelude::*;
use std::fs::File;
use ::plot::*;

const CHR_WIDTH: f64 = 4.0;

pub struct FlatPlotter {
    result: RunResult,
    settings: Settings,

    max_length: f64,
    width: f64,
    height: f64,
}

impl Plotter for FlatPlotter {
    fn new(settings: Settings, result: RunResult) -> FlatPlotter {
        let length = result.strand.length as f64;
        FlatPlotter {
            result,
            settings,

            max_length: length,
            width: 1500.0,
            height: 230.0,
        }
    }

    fn plot(&self) -> Result<()> {
        let out_filename = format!("{}.svg", &self.settings.out_file);
        File::create(&out_filename)
            .and_then(|mut f| f.write_all(self.plot_flat().as_bytes()))
            .and_then(|_| { println!("Flat plot written to `{}`", &out_filename); Ok(()) })
            .chain_err(|| format!("Unable to write in `{}`", &out_filename))?;

        Ok(())
    }
}

impl FlatPlotter {
    fn plot_flat(&self) -> String {
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
                        self.result.strand.length as f64/self.max_length*self.width,
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
                        self.result.strand.length as f64/self.max_length*self.width,
                        self.height-CHR_WIDTH/2.0,
                        CHR_WIDTH
        );

        //
        // Ticks
        //
        for i in 0..self.max_length as i64 {
            if i % 1_000_000 == 0 {
                let height = if i % 10_000_000 == 0 {
                    self.height + 15.0
                } else if i % 5_000_000 == 0 {
                    self.height + 10.0
                } else {
                    self.height + 5.0
                };
                let x = i as f64/self.max_length*self.width;

                svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#898989' stroke-width='1'/>",
                                x, self.height,
                                x, height
                );

                if i % 10_000_000 == 0 {
                    svg += &format!("<text x='{}' y='{}' font-family='Helvetica' font-size='8'>{}Mb</text>",
                                    x+4.0,
                                    height,
                                    i / 1_000_000);
                }
            }
        }

        //
        // Feaures
        //
        for features_family in &self.settings.feature_tracks {
            for feature in features_family.iter() {
                for position in &feature.positions {
                    let (start, end) = match *position {
                        FeaturePosition::Relative { ref chr, start, length} => {
                            let chr = self.result.strand.find_chr(&chr).unwrap_or_else(|| panic!("Unable to find fragment `{}`", chr));
                            (chr.position + start, chr.position + start + length)
                        }
                        FeaturePosition::Absolute { start, length }         => { (start, start + length) }
                    };

                    let color = format!("#{:2X}{:2X}{:2X}", rand::random::<i8>(), rand::random::<i8>(), rand::random::<i8>());
                    let x0 = start as f64/self.max_length * self.width;
                    let x1 = end as f64/self.max_length * self.width;
                    let x2 = x1 + 2.0;
                    let x3 = x0 - 2.0;
                    let font_size = 8.0;

                    svg += &format!("<polygon points='{},{} {},{} {},{} {},{}' style='fill:{};'/>\n",
                                    x0, self.height,
                                    x1, self.height,
                                    x2, self.height + 10.0,
                                    x3, self.height + 10.0,
                                    color
                    );

                    svg += &format!("<text x='{}' y='{}' font-family='sans-serif' font-size='{}' style='writing-mode: tb;'>{}</text>",
                                    x0, self.height + 20.0 + font_size,
                                    font_size, feature.name);
                }
            }
        }


        for family in &self.result.families {
            for sd in family {
                let left1 = (sd.global_left_position as f64)/self.max_length * self.width;
                let left2 = (sd.global_left_position as f64 + sd.left_length as f64)/self.max_length * self.width;
                let right1 = (sd.global_right_position as f64)/self.max_length * self.width;
                let right2 = (sd.global_right_position as f64 + sd.right_length as f64)/self.max_length * self.width;

                let color = if sd.reversed { &self.settings.color2 } else { &self.settings.color1 };

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
                                if left2 - left1 < self.settings.min_thickness { left1 + self.settings.min_thickness} else { left2 }, CHR_WIDTH,
                                if right2 - right1 < self.settings.min_thickness { right1 + self.settings.min_thickness} else { right2 }, self.height - CHR_WIDTH,
                                right1, self.height - CHR_WIDTH,
                                color,
                                color,
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

        format!("<?xml version='1.0' encoding='UTF-8' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>{}</svg>",
                self.width,
                self.height+300.0,
                svg)
    }
}
