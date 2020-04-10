use crate::plot::*;
use crate::plot::colorizers::Colorizer;

use std::fs::File;
use std::io::Write;

use thousands::Separable;

pub struct GenomePlotter {
    result: RunResult,
    settings: Settings,
    colorizer: Box<dyn Colorizer>,
}

impl Plotter for GenomePlotter {
    fn new(settings: Settings, result: RunResult, colorizer: Box<dyn Colorizer>) -> GenomePlotter {
        GenomePlotter {
            result,
            settings,
            colorizer,
        }
    }

    fn plot(&self) -> Result<()> {
        let out_filename = format!("{}.svg", &self.settings.out_file);
        File::create(&out_filename)
            .and_then(|mut f| f.write_all(self.plot_genome().as_bytes()))
            .and_then(|_| { log::info!("Genome plot written to `{}`", &out_filename); Ok(()) })
            .with_context(|| format!("Failed to save plot to `{}`", &out_filename))?;

        Ok(())
    }
}

impl GenomePlotter {
    fn plot_genome(&self) -> String {
        let mut svg = String::new();

        // 100 px/chromosome
        let chr_spacing = 100.0;
        let chr_width = 40.0;
        let height_factor = 800.0;
        let factor = 1.0/self.result.strand.map.iter().map(|chr| chr.length).max().unwrap() as f64 * height_factor;

        let width = chr_spacing as usize * (self.result.strand.map.len() + 1);
        // height_factor + 50px top + 100px bot
        let height = height_factor + 50.0 + 100.0;

        // 1. Draw the chromosomes
        for (i, chr) in self.result.strand.map.iter().enumerate() {
            // Chromosome bar
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}44' stroke-width='{}'/>\n",
                            chr_spacing + i as f64 * chr_spacing,
                            50,
                            chr_spacing + i as f64 * chr_spacing,
                            50.0 + factor*chr.length as f64,
                            self.colorizer.color_fragment(&chr.name),
                            chr_width,
            );

            // Central delimiter
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#111' stroke-width='{}' stroke-dasharray='5,5'/>\n",
                            chr_spacing + i as f64 * chr_spacing,
                            50,
                            chr_spacing + i as f64 * chr_spacing,
                            50.0 + factor*chr.length as f64,
                            1,
            );
            // Side delimiters
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#222' stroke-width='{}' stroke-dasharray='1,2'/>\n",
                            chr_spacing + i as f64 * chr_spacing - chr_width/4.0,
                            50,
                            chr_spacing + i as f64 * chr_spacing - chr_width/4.0,
                            50.0 + factor*chr.length as f64,
                            0.5,
            );
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#222' stroke-width='{}' stroke-dasharray='1,2'/>\n",
                            chr_spacing + i as f64 * chr_spacing + chr_width/4.0,
                            50,
                            chr_spacing + i as f64 * chr_spacing + chr_width/4.0,
                            50.0 + factor*chr.length as f64,
                            0.5,
            );

            // Label
            svg += &format!("<text x='{}' y='{}' style='font-size: 11;'>{}</text>\n",
                            chr_spacing + i as f64 * chr_spacing - 10.0,
                            20 + (i%2) * 10,
                            if chr.name.len() > 8 { &chr.name[0..3] } else { &chr.name },
            );
        }

        // 2. Draw the SDs
        for family in &self.result.families {
            for sd in family {
                let color = self.colorizer.color(sd);
                let x: Box<dyn Fn(usize) -> f64> = match (sd.chr_left == sd.chr_right, sd.reversed) {
                    (true,  false) => {Box::new(|x| chr_spacing - 3.0*chr_width/8.0 + chr_spacing*x as f64)}
                    (true,  true)  => {Box::new(|x| chr_spacing - 1.0*chr_width/8.0 + chr_spacing*x as f64)}
                    (false, false) => {Box::new(|x| chr_spacing + 1.0*chr_width/8.0 + chr_spacing*x as f64)}
                    (false, true)  => {Box::new(|x| chr_spacing + 3.0*chr_width/8.0 + chr_spacing*x as f64)}
                };

                // left arm
                if sd.chr_left != COLLAPSED_NAME {
                    let chr_left_index = self.result.strand.find_chr_index(&sd.chr_left).unwrap();
                    let left = sd.chr_left_position;
                    let start = factor * left as f64;
                    let mut end = factor * (left + sd.left_length) as f64;
                    if start - end < self.settings.min_thickness {
                        end = start + self.settings.min_thickness
                    };
                    svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='{}'><title>{}</title></line>\n",
                                    x(chr_left_index),
                                    50.0 + start,
                                    x(chr_left_index),
                                    50.0 + end,
                                    color,
                                    chr_width/4.0,
                                    &format!("{}: {} → {}  ({}bp)\n{}: {} → {} ({}bp)",
                                             &sd.chr_left,
                                             sd.chr_left_position.separate_with_spaces(),
                                             (sd.chr_left_position+sd.left_length).separate_with_spaces(),
                                             sd.left_length.separate_with_spaces(),

                                             &sd.chr_right,
                                             sd.chr_right_position.separate_with_spaces(),
                                             (sd.chr_right_position+sd.right_length).separate_with_spaces(),
                                             sd.right_length.separate_with_spaces()
                                    )
                    );
                }

                // right arm
                if sd.chr_right != COLLAPSED_NAME {
                    let chr_right_index = self.result.strand.find_chr_index(&sd.chr_right).unwrap();
                    let right = sd.chr_right_position;
                    let start = factor * right as f64;
                    let mut end = factor * (right + sd.right_length) as f64;
                    if start - end < self.settings.min_thickness {
                        end = start + self.settings.min_thickness
                    };
                    svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='{}'><title>{}</title></line>\n",
                                    x(chr_right_index),
                                    50.0 + start,
                                    x(chr_right_index),
                                    50.0 + end,
                                    color,
                                    chr_width/4.0,
                                    &format!("{}: {} → {}  ({}bp)\n{}: {} → {} ({}bp)",
                                             &sd.chr_left,
                                             sd.chr_left_position.separate_with_spaces(),
                                             (sd.chr_left_position+sd.left_length).separate_with_spaces(),
                                             sd.left_length.separate_with_spaces(),

                                             &sd.chr_right,
                                             sd.chr_right_position.separate_with_spaces(),
                                             (sd.chr_right_position+sd.right_length).separate_with_spaces(),
                                             sd.right_length.separate_with_spaces()
                                    )

                    );
                }
            }
        }

        // 3. Label everything

        // 4. Done!
        format!(
            r#"
<!DOCTYPE svg PUBLIC '-//W3C//DTD SVG 1.0//EN' 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'>
<svg version='1.0' width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>
{}
</svg>"#,
            width,
            height,
            svg)
    }
}
