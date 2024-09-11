use std::fs::File;
use std::io::prelude::*;
use thousands::Separable;

use crate::plot::colorizers::Colorizer;
use crate::plot::*;

#[derive(Clone)]
enum SpanClass {
    Duplicon {
        reversed: bool,
        complemented: bool,
        both: bool,
    },
    Feature,
}
#[derive(Clone)]
struct Span {
    start: usize,
    length: usize,
    class: SpanClass,
}

#[derive(Debug)]
enum DrawCommand {
    Distance(i64),
    Feature { length: i64, color: String },
}

pub struct RosaryPlotter {
    result: RunResult,
    settings: Settings,
    colorizer: Box<dyn Colorizer>,

    clustering_margin: usize,
    rosary_mode: bool,
}

impl Plotter for RosaryPlotter {
    fn plot(&self) -> Result<()> {
        let out_filename = format!("{}.svg", &self.settings.out_file);
        File::create(&out_filename)
            .and_then(|mut f| f.write_all(self.plot_squish().as_bytes()))
            .map(|_| {
                log::info!("Rosary plot written to `{}`", &out_filename);
                
            })
            .with_context(|| format!("Failed to save plot to `{}`", &out_filename))?;

        Ok(())
    }
}

impl RosaryPlotter {
    pub fn new(
        settings: Settings,
        result: RunResult,
        colorizer: Box<dyn Colorizer>,
        clustering_margin: usize,
        rosary_mode: bool,
    ) -> RosaryPlotter {
        log::info!("Clustering margin: {}bp", clustering_margin);
        RosaryPlotter {
            result,
            settings,
            colorizer,

            clustering_margin,
            rosary_mode,
        }
    }

    fn size_for_feature(l: f32) -> f32 {
        l / 10_000.
    }

    fn size_for_void(l: f32) -> f32 {
        (l / 100_000.).sqrt()
    }

    fn annotations_for_chr(&self, chr: &Start) -> Vec<Span> {
        self.settings
            .feature_tracks
            .iter()
            .flat_map(|family| {
                family.iter().flat_map(|feature| {
                    feature.positions.iter().filter_map(|position| {
                        match position {
                            FeaturePosition::Relative {
                                chr: my_chr,
                                start,
                                length,
                            } => {
                                let my_chr =
                                    self.result.strand.find_chr(my_chr).unwrap_or_else(|| {
                                        panic!("Unable to find fragment `{}`", my_chr)
                                    });
                                if my_chr.name == chr.name {
                                    Some(Span {
                                        start: *start,
                                        length: *length,
                                        class: SpanClass::Feature,
                                    })
                                } else {
                                    None
                                }
                            }
                            FeaturePosition::Absolute {
                                start: _,
                                length: _,
                            } =>
                            // (*start, start + length),
                            {
                                unimplemented!()
                            }
                        }
                    })
                })
            })
            .collect()
    }

    fn duplicons_for_chr(&self, chr: &Start) -> Vec<Span> {
        let mut proto_duplicons: Vec<Span> = self
            .result
            .families
            .iter()
            .flatten()
            .filter(|sd| sd.chr_left == chr.name || sd.chr_right == chr.name)
            .flat_map(|sd| {
                vec![
                    (
                        sd.chr_left.to_owned(),
                        sd.chr_left_position,
                        sd.left_length,
                        sd.reversed,
                        sd.complemented,
                    ),
                    (
                        sd.chr_right.to_owned(),
                        sd.chr_right_position,
                        sd.right_length,
                        sd.reversed,
                        sd.complemented,
                    ),
                ]
            })
            .filter(|duplicon| duplicon.0 == chr.name)
            .map(|duplicon| Span {
                start: duplicon.1,
                length: duplicon.2,
                class: SpanClass::Duplicon {
                    reversed: duplicon.3,
                    complemented: duplicon.4,
                    both: false,
                },
            })
            .collect();
        proto_duplicons.sort_by_key(|span| span.start);

        let mut duplicons: Vec<Span> = Vec::new();
        for new in proto_duplicons {
            if let Some(last) = duplicons.last_mut() {
                if new.start <= last.start + last.length + self.clustering_margin {
                    last.length = new.start + new.length - last.start;
                    if let SpanClass::Duplicon {
                        reversed: old_r,
                        complemented: old_c,
                        ref mut both,
                    } = last.class
                    {
                        if let SpanClass::Duplicon {
                            reversed: new_r,
                            complemented: new_c,
                            ..
                        } = new.class
                        {
                            if (old_r != new_r) || (old_c != new_c) {
                                *both = true;
                            }
                        }
                    }
                } else {
                    duplicons.push(new);
                }
            } else {
                duplicons.push(new);
            }
        }

        duplicons
    }

    fn plot_squish(&self) -> String {
        let chr_draw_commands: Vec<Vec<DrawCommand>> = self
            .result
            .strand
            .map
            .iter()
            .map(|current_chr| {
                let duplicons = self.duplicons_for_chr(current_chr);
                let annotations = self.annotations_for_chr(current_chr);

                log::trace!(
                    "{} :: Plotting {} duplicons and {} features",
                    current_chr.name,
                    duplicons.len(),
                    annotations.len()
                );
                let mut features = duplicons
                    .into_iter()
                    .chain(annotations)
                    .collect::<Vec<Span>>();
                features.sort_by_key(|f| f.start);

                let mut pos = 0;
                let mut draw_commands = Vec::<DrawCommand>::new();
                for span in features {
                    let mut distance = span.start - pos;
                    if self.rosary_mode {
                        while distance > 0 {
                            if distance > 10_000_000 {
                                draw_commands.push(DrawCommand::Distance(10_000_000));
                                distance -= 10_000_000;
                            } else if distance > 1_000_000 {
                                draw_commands.push(DrawCommand::Distance(1_000_000));
                                distance -= 1_000_000;
                            } else if distance > 100_000 {
                                draw_commands.push(DrawCommand::Distance(100_000));
                                distance -= 100_000;
                            } else {
                                draw_commands.push(DrawCommand::Distance(distance as i64));
                                distance = 0;
                            }
                        }
                    } else {
                        draw_commands.push(DrawCommand::Distance(distance as i64));
                    }

                    draw_commands.push(DrawCommand::Feature {
                        length: span.length as i64,
                        color: if let SpanClass::Duplicon {
                            reversed,
                            complemented,
                            both,
                            ..
                        } = span.class
                        {
                            if !both {
                                if reversed && complemented {
                                    "#00b2ae".to_owned()
                                } else {
                                    "#ff5b00".to_owned()
                                }
                            } else {
                                "#9741ad".to_owned()
                            }
                        } else {
                            "#66491e".to_owned()
                        },
                    });
                    pos = span.start + span.length;
                }
                if pos < current_chr.length {
                    draw_commands.push(DrawCommand::Distance((current_chr.length - pos) as i64));
                }
                draw_commands
            })
            .collect();

        const SCALES: &[(i64, &str)] = &[
            (100_000, "100kbp"),
            (1_000_000, "1Mbp"),
            (5_000_000, "5Mbp"),
            (10_000_000, "10Mbp"),
            (50_000_000, "50Mbp"),
        ];
        let largest_bead = chr_draw_commands
            .iter()
            .flat_map(|x| x.iter())
            .filter_map(|c| match c {
                DrawCommand::Feature { .. } => None,
                DrawCommand::Distance(l) => Some(*l),
            })
            .max()
            .unwrap_or(0);

        let largest_square = chr_draw_commands
            .iter()
            .flat_map(|x| x.iter())
            .filter_map(|c| match c {
                DrawCommand::Feature { length: l, .. } => Some(*l),
                DrawCommand::Distance(..) => None,
            })
            .max()
            .unwrap_or(0);

        let captions_beads_text = SvgObject::Text {
            x: 0.,
            y: 0.,
            text: "Duplications-devoid regions".to_string(),
            font_size: None,
            color: None,
        };
        let captions_beads = SCALES
            .iter()
            .filter(|x| x.0 <= largest_bead)
            .fold(
                (
                    (0., captions_beads_text.dims().1 + 5.),
                    SvgGroup::new().push(captions_beads_text),
                ),
                |((x, y), g), l| {
                    let r = Self::size_for_void(l.0 as f32);
                    let text = SvgObject::Text {
                        x,
                        y,
                        text: l.1.to_string(),
                        color: None,
                        font_size: None,
                    };
                    let bead = SvgObject::Circle {
                        cx: x + text.dims().0 / 3.,
                        cy: y + text.dims().1 + 5.,
                        r,
                        fill: "#555555".to_string(),
                    };
                    (
                        (x + text.dims().0 + bead.dims().0 + 10., y),
                        g.append(SvgGroup::new().push(bead).push(text)),
                    )
                },
            )
            .1;

        let captions_squares_text = SvgObject::Text {
            x: 0.,
            y: 0.,
            text: "Duplications-rich regions".to_string(),
            font_size: None,
            color: None,
        };
        let captions_squares = SCALES
            .iter()
            .filter(|x| x.0 <= largest_square)
            .fold(
                (
                    (0., captions_squares_text.dims().1 + 5.),
                    SvgGroup::new().push(captions_squares_text),
                ),
                |((x, y), g), l| {
                    let w = Self::size_for_feature(l.0 as f32);
                    let text = SvgObject::Text {
                        x,
                        y,
                        text: l.1.to_string(),
                        color: None,
                        font_size: None,
                    };

                    let square = SvgObject::Line {
                        x1: x + text.dims().0 / 3.,
                        x2: x + text.dims().0 / 3.,
                        y1: y + text.dims().1 + 5.,
                        y2: y + text.dims().1 + w + 5.,
                        stroke: Some("#bbb".to_string()),
                        stroke_width: w,
                        hover: None,
                    };
                    (
                        (x + text.dims().0 + square.dims().0 + 10., y),
                        g.append(SvgGroup::new().push(square).push(text)),
                    )
                },
            )
            .1;

        let captions = SvgGroup::new()
            .append(captions_squares.shift(0., captions_beads.dims().1 + 15.))
            .append(captions_beads);

        let labels = self
            .result
            .strand
            .map
            .iter()
            .map(|chr| SvgObject::Text {
                x: 0.,
                y: 0.,
                text: chr.name.clone(),
                color: None,
                font_size: None,
            })
            .collect::<Vec<_>>();

        let label_space = 5.0
            + labels
                .iter()
                .map(|x| (x.dims().0 + 1.) as i64)
                .max()
                .unwrap() as f32;

        let chrs: Vec<SvgGroup> = chr_draw_commands
            .iter()
            .map(|chr_commands| {
                chr_commands
                    .iter()
                    .fold((label_space, SvgGroup::new()), |(x, ax), cmd| match cmd {
                        DrawCommand::Distance(d) => {
                            let r = Self::size_for_void(*d as f32);
                            let cx = x + r;
                            let cy = 0.;
                            (
                                x + 2. * r,
                                ax.push(SvgObject::Circle {
                                    cx,
                                    cy,
                                    r,
                                    fill: "#555555".to_string(),
                                }),
                            )
                        }
                        DrawCommand::Feature { length, color } => {
                            let width = Self::size_for_feature(*length as f32);
                            let x1 = x;
                            let y1 = 0.;
                            let x2 = x1 + width;
                            let y2 = 0.;
                            (
                                x2,
                                ax.push(SvgObject::Line {
                                    x1,
                                    y1,
                                    x2,
                                    y2,
                                    stroke: Some(color.to_string()),
                                    stroke_width: width,
                                    hover: Some(format!(
                                        "{} → {}  ({}bp)",
                                        "na",
                                        "na",
                                        length.separate_with_spaces()
                                    )),
                                }),
                            )
                        }
                    })
                    .1
            })
            .collect();

        let main_plot = labels
            .into_iter()
            .zip(chrs)
            .fold((0., SvgGroup::new()), |(y, ax), (mut label, chr)| {
                let height = label.dims().1.max(chr.dims().1);
                let shift = y + height / 2.0;
                label.shift(0., shift);
                (
                    y + height + 10.,
                    ax.push(label).append(chr.shift(0., shift)),
                )
            })
            .1
            .shift(0., captions.dims().1 + 20.);
        let all = SvgGroup::new()
            .append(captions)
            .append(main_plot)
            .shift(10., 10.)
        // .transpose()
            ;

        format!(
            "<?xml version='1.0' encoding='UTF-8'  standalone='no' ?> <!DOCTYPE svg \
             PUBLIC '-//W3C//DTD SVG 1.0//EN' 'http://www.w3.org/TR/2001/REC-SVG-20010904/DTD/svg10.dtd'> \
             <svg version='1.0' width='{}' height='{}' xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink'>\n \
             {} \
             </svg>",
            all.dims().0 + 10., all.dims().1 + 10.,
            all.render(),
        )
    }
}
