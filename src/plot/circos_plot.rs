use crate::plot::*;
use crate::plot::colorizers::Colorizer;
use crate::utils::slugify;

use std::fs::File;
use std::io::prelude::*;

pub struct CircosPlotter {
    result: RunResult,
    settings: Settings,
}

impl Plotter for CircosPlotter {
    fn new(settings: Settings, result: RunResult, _x: Box<dyn Colorizer>) -> CircosPlotter {
        CircosPlotter {
            result,
            settings,
        }
    }

    fn plot(&self) -> Result<()> {
        let prefix = self.settings.out_file.clone();

        let karyotype_filename = &format!("{}.karyotype", &prefix);
        let links_filename = &format!("{}.links", &prefix);
        let config_filename = &format!("{}.conf", &prefix);

        File::create(karyotype_filename)
            .and_then(|mut f| f.write_all(self.plot_karyotype().as_bytes()))
            .and_then(|_| {
                log::info!("Karyotype written to `{}`", &karyotype_filename);
                Ok(())
            })
            .with_context(|| format!("Failed to save karyotype to `{}`", &karyotype_filename))?;

        File::create(links_filename)
            .and_then(|mut f| f.write_all(self.plot_links().as_bytes()))
            .and_then(|_| {
                log::info!("Links written to `{}`", &links_filename);
                Ok(())
            })
            .with_context(|| format!("Failed to save links to `{}`", &links_filename))?;

        File::create(config_filename)
            .and_then(|mut f| {
                f.write_all(
                    self.plot_config(karyotype_filename, links_filename)
                        .as_bytes(),
                )
            })
            .and_then(|_| {
                log::info!("Config written to `{}`", &config_filename);
                Ok(())
            })
            .with_context(|| format!("Failed to save plot configuration to `{}`", &config_filename))?;

        println!(
            "\nYou can now edit `{}` and/or run `circos {}` to generate the final plot.",
            config_filename, config_filename
        );
        Ok(())
    }
}

impl CircosPlotter {
    fn plot_karyotype(&self) -> String {
        fn encode_chromosome(chr: &Start) -> String {
            format!(
                "chr - {id} {label} {start} {end} {color}",
                id = slugify(&chr.name),
                label = slugify(&chr.name),
                start = 0,
                end = chr.length,
                color = "grey"
            )
        }

        self.result
            .strand
            .map
            .iter()
            .map(encode_chromosome)
            .collect::<Vec<String>>()
            .join("\n")
    }

    fn plot_links(&self) -> String {
        let links : String = self.result.families
            .iter()
            .map(|family|
                 family
                 .iter()
                 .map(
                     |sd|
                     format!("{chr_left} {chr_left_start} {chr_left_end} {chr_right} {chr_right_start} {chr_right_end} {color}",
                             chr_left        = slugify(&sd.chr_left),
                             chr_left_start  = sd.chr_left_position,
                             chr_left_end    = sd.chr_left_position + sd.left_length,
                             chr_right       = slugify(&sd.chr_right),
                             chr_right_start = sd.chr_right_position,
                             chr_right_end   = sd.chr_right_position + sd.right_length,
                             color           = if sd.reversed { "color=teal" } else { "color=orange" }
                     )
                 )
            )
            .flatten()
            .collect::<Vec<String>>()
            .join("\n");

        links
    }

    fn plot_config(&self, karyotype_filename: &str, links_filename: &str) -> String {
        format!("
karyotype = {karyotype_filename}
chromosomes_units = 1000000

<colors>
orange = 255,  91,   0, 0.5
teal   =   0, 178, 174, 0.5
</colors>

### IDEOGRAM SECTION
<ideogram>

<spacing>
default = 0.005r
</spacing>

radius           = 0.90r
thickness        = 20p
fill             = yes
stroke_color     = dgrey
stroke_thickness = 2p
show_label       = yes
label_font       = default
label_radius     = dims(image,radius) - 60p
label_size       = 30
label_parallel   = yes

</ideogram>
### END IDEOGRAM SECTION

### TICKS SECTION
show_ticks          = yes
show_tick_labels    = yes

<ticks>
radius           = 1r
color            = black
thickness        = 2p
multiplier       = 1e-6
format           = %d

<tick>
spacing        = 5u
size           = 10p
</tick>

<tick>
spacing        = 25u
size           = 15p
show_label     = yes
label_size     = 20p
label_offset   = 10p
format         = %d
</tick>
</ticks>
### END TICKS SECTION

<links>
   <link>
      file          = {links_filename}
      radius        = 0.95r
      bezier_radius = 0r
      ribbon        = yes
   </link>
</links>

<image>
<<include {circos_root}/etc/image.conf>>
</image>
<<include {circos_root}/etc/colors_fonts_patterns.conf>>
<<include {circos_root}/etc/housekeeping.conf>>
",
                karyotype_filename = karyotype_filename,
                links_filename = links_filename,
                circos_root = std::env::var("CIRCOS_ROOT").unwrap_or_else(|_| {
                    log::error!("CIRCOS_ROOT is not set - using a placeholder in config file.\nPlease set the CICRCOS_ROOT environment variable with the root of your Circos install.");
                    "REPLACE_ME_WITH_CIRCOS_ROOT".to_string()
                }),
        )
    }
}
