extern crate rand;
use plot::*;
use plot::colorizers::Colorizer;

use std::io::prelude::*;
use std::fs::File;
use separator::Separatable;

const CHR_WIDTH: f64 = 4.0;

pub struct FlatPlotter {
    result: RunResult,
    settings: Settings,
    colorizer: Box<dyn Colorizer>,

    max_length: f64,
    width: f64,
    height: f64,
}

impl Plotter for FlatPlotter {
    fn new(settings: Settings, result: RunResult, colorizer: Box<dyn Colorizer>) -> FlatPlotter {
        let length = result.strand.length as f64;
        FlatPlotter {
            result,
            settings,
            colorizer,

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

        let mut offset: i64 = 0;
        for (j, chr) in self.result.strand.map.iter().enumerate() {
            //
            // Chromosomes
            //
            // Top bar
            svg += &format!(
                "<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='{}'/>",
                offset as f64/self.max_length*self.width, CHR_WIDTH/2.0,
                (offset + chr.length as i64) as f64/self.max_length*self.width,
                CHR_WIDTH/2.0,
                self.colorizer.color_fragment(&chr.name), CHR_WIDTH
            );
            // Bottom bar
            svg += &format!(
                "<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='{}'/>",
                offset as f64/self.max_length*self.width, self.height-CHR_WIDTH/2.0,
                (offset + chr.length as i64) as f64/self.max_length*self.width,
                self.height-CHR_WIDTH/2.0,
                self.colorizer.color_fragment(&chr.name), CHR_WIDTH
            );
            // Name
            svg += &format!(
                "<text x='{}' y='{}' font-family='Helvetica' font-size='12'>{}</text>",
                offset as f64/self.max_length*self.width, self.height + 35.0, chr.name
            );

            //
            // Ticks
            //
            for i in 0..chr.length as i64 {
                if i % 1_000_000 == 0 {
                    let height = if i % 10_000_000 == 0 {
                        self.height + 7.0
                    } else if i % 5_000_000 == 0 {
                        self.height + 5.0
                    } else {
                        self.height + 3.0
                    };
                    let x = (i + offset) as f64/self.max_length*self.width;

                    svg += &format!(
                        "<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#898989' stroke-width='1'/>",
                        x, self.height,
                        x, height
                    );

                    if i % 10_000_000 == 0 {
                        svg += &format!(
                            "<text x='{}' y='{}' font-family='Helvetica' font-size='8'>{}Mb</text>",
                            x,
                            self.height + 15.0 + ((j % 2) as f64)*5.0,
                            i / 1_000_000
                        );
                    }
                }
            }

            offset += chr.length as i64;
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

                let color = self.colorizer.color(sd);

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
                                &format!("{}: {} → {}  ({}bp)\n{}: {} → {} ({}bp)",
                                         &sd.chr_left,
                                         sd.chr_left_position.separated_string(),
                                         (sd.chr_left_position+sd.left_length).separated_string(),
                                         sd.left_length.separated_string(),

                                         &sd.chr_right,
                                         sd.chr_right_position.separated_string(),
                                         (sd.chr_right_position+sd.right_length).separated_string(),
                                         sd.right_length.separated_string()
                                )
                );
            }
        }

        format!("<?xml version='1.0' encoding='UTF-8' standalone='no' ?> <!DOCTYPE svg \
                 PUBLIC '-//W3C//DTD SVG 1.0//EN' \
                 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> <svg version='1.0' \
                 width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' \
                 xmlns:xlink='http://www.w3.org/1999/xlink'>{}</svg>",
                self.width+25.0,
                self.height+40.0,
                svg)
    }
}
