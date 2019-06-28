use ::plot::*;

use std::fs::File;
use std::io::Write;

pub struct GenomePlotter {
    result: RunResult,
    settings: Settings,
}

impl Plotter for GenomePlotter {
    fn new(settings: Settings, result: RunResult) -> GenomePlotter {
        GenomePlotter {
            result,
            settings,
        }
    }

    fn plot(&self) -> Result<()> {
        let out_filename = format!("{}.svg", &self.settings.out_file);
        File::create(&out_filename)
            .and_then(|mut f| f.write_all(self.plot_genome().as_bytes()))
            .and_then(|_| { println!("Genome plot written to `{}`", &out_filename); Ok(()) })
            .chain_err(|| format!("Unable to write in `{}`", &out_filename))?;

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
        let factor = 1.0/self.result.strand1.map.iter().map(|chr| chr.length).max().unwrap() as f32 * height_factor;
        let threshold = 0.1;

        let width = chr_spacing as usize * (self.result.strand1.map.len() + 1);
        // height_factor + 50px top + 100px bot
        let height = height_factor + 50.0 + 100.0;


        // 1. Draw the chromosomes
        // let regex_chr_name = Regex::new("(?:chromosome ?([0-9XYxyZWzw]{1,2}))|(?:chr ?([0-9XYxyZWzwA-Za-z]{1,2}))|(?:[^0-9XYxyZWzwA-Za-z]([0-9XYxyZWzw]{1,2})[^0-9XYxyZWzw.:|,a-zA-Z])")
        //     .expect("Unable to create regex");
        for (i, chr) in self.result.strand1.map.iter().enumerate() {
            // Chromosome bar
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#ccc' stroke-width='{}'/>\n",
                            chr_spacing + i as f32 * chr_spacing,
                            50,
                            chr_spacing + i as f32 * chr_spacing,
                            50.0 + factor*chr.length as f32,
                            chr_width,
            );

            // Central delimiter
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#111' stroke-width='{}' stroke-dasharray='5,5'/>\n",
                            chr_spacing + i as f32 * chr_spacing,
                            50,
                            chr_spacing + i as f32 * chr_spacing,
                            50.0 + factor*chr.length as f32,
                            1,
            );
            // Side delimiters
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#222' stroke-width='{}' stroke-dasharray='1,2'/>\n",
                            chr_spacing + i as f32 * chr_spacing - chr_width/4.0,
                            50,
                            chr_spacing + i as f32 * chr_spacing - chr_width/4.0,
                            50.0 + factor*chr.length as f32,
                            0.5,
            );
            svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#222' stroke-width='{}' stroke-dasharray='1,2'/>\n",
                            chr_spacing + i as f32 * chr_spacing + chr_width/4.0,
                            50,
                            chr_spacing + i as f32 * chr_spacing + chr_width/4.0,
                            50.0 + factor*chr.length as f32,
                            0.5,
            );


            // Label
            // let chr_name_cap = regex_chr_name.captures(&chr.name).expect("No matches");
            // let chr_name = chr_name_cap.get(1).unwrap().as_str();
            svg += &format!("<text x='{}' y='{}' style='font-size: 11;'>{}</text>\n",
                            chr_spacing + i as f32 * chr_spacing - 10.0,
                            20 + (i%2) * 10,
                            if chr.name.len() > 8 { &chr.name[0..3] } else { &chr.name },
            );
        }

        // 2. Draw the SDs
        for family in &self.result.families {
            for sd in family {
                let color = if sd.reversed { &self.settings.color2 } else  { &self.settings.color1 };
                let x: Box<dyn Fn(usize) -> f32> = match (sd.chr_left == sd.chr_right, sd.reversed) {
                    (true,  false) => {Box::new(|x| chr_spacing - 3.0*chr_width/8.0 + chr_spacing*x as f32)}
                    (true,  true)  => {Box::new(|x| chr_spacing - 1.0*chr_width/8.0 + chr_spacing*x as f32)}
                    (false, false) => {Box::new(|x| chr_spacing + 1.0*chr_width/8.0 + chr_spacing*x as f32)}
                    (false, true)  => {Box::new(|x| chr_spacing + 3.0*chr_width/8.0 + chr_spacing*x as f32)}
                };

                // left arm
                let chr_left_index = self.result.strand1.find_chr_index(&sd.chr_left).unwrap();
                let left = sd.chr_left_position;
                let start = factor * left as f32;
                let mut end = factor * (left + sd.length) as f32;
                if start - end < threshold { end += threshold };
                svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='{}'/>\n",
                                x(chr_left_index),
                                50.0 + start,
                                x(chr_left_index),
                                50.0 + end,
                                color,
                                chr_width/4.0,
                );

                // right arm
                let chr_right_index = self.result.strand2.find_chr_index(&sd.chr_right).unwrap();
                let right = sd.chr_right_position;
                let start = factor * right as f32;
                let mut end = factor * (right + sd.length) as f32;
                if start - end < threshold { end += threshold };
                svg += &format!("<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='{}'/>\n",
                                x(chr_right_index),
                                50.0 + start,
                                x(chr_right_index),
                                50.0 + end,
                                color,
                                chr_width/4.0,
                );
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
