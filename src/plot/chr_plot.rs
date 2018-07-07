use ::structs::*;
use ::plot::{Plotter, Settings}

pub struct GenomePlotter {
    result: RunResult,
    settings: Settings,
}

impl Plotter for GenomePlotter {
    fn new(settings: Settings, result: RunResult) -> GenomePlotter {
        GenomePlotter {
            result: result,
            settings: settings,
        }
    }

    fn plot(self) {
        let mut f = File::create(&self.setings.out_file).unwrap();
        f.write_all(self.plot().as_bytes()).expect("Unable to write result to file");
    }
}

impl GenomPlotter {
    fn plot(self) -> String {
        let mut svg = String::new();

        // 100 px/chromosome
        let width = 100 * self.result.map.length();
        // 500px + 50px top + 100px bot
        let height = 650;
        let factor = 1.0/ * 500.0;


        // 1. Draw the chromosomes
        for (i, chr) in self.result.map.iter().enumerate() {
            svg += &format!(r#"<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='#ccc' stroke-width='20'/>"#,
                            50 + i*100,
                            50,
                            50 + i*100,
                            factor*chr.length
            )
        }

        // 2. Draw the SDs
        for &sd in self.result.sds {
            let (chr_left, chr_right) = sd.chr_index();
            let x = if inter { |x| { 45 + 100*x }} else { |x| { 55 + 100*x } }
            // left arm
            svg += &format!(r#"<line x1='{}' y1='{}' x2='{}' y2='{}' stroke='{}' stroke-width='10'"#,
                            x(chr_left),
                            50,
                            x(chr_left),
                            factor * sd.left/chr_left

            );
            // right arm
        }

        // 3. Label eveyrthing

        // 4. Done!
    }
}
